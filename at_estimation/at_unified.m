function at = at_unified(X, varargin)

% Estimate activation times (AT) from transmembrane voltages (TMV) or 
% extracellular potentials (EP).
% Depending on the input, AT estimation is based on detecting the
% intrinsic deflection of a derivative signal (power == inf)
% or on fitting ATs to inter-node delays obtained using
% cross-correlation of the derivative signal (power < inf).
%
% at = at_unified(X, NAME-VALUE PAIRS)
%
% INPUTS:
%              X: TMVs or EPs [numNodes x numTimesteps].
%                 The positive slope is detected, so EPs must be inverted.
%
%                 NAME-VALUE PAIRS:
%
%   'upsampling': Factor for temporal upsampling (linear interpolation).
%                 Upsampling is particularly important to capture small
%                 time shifts between signals in integer-type delay values
%                 obtained using cross-correlation.
%                 Default: 10
%
%   'derivative': Type of derivative to be used:
%                 't':  temporal,
%                 's':  spatial,
%                 'st': spatiotemporal.
%                 Default: 't'
%
%        'sigma': Std. dev. in samples used for temporal Gaussian filtering
%                 (zero-phase, 2nd order).
%                 Temporal filtering only affects the temporal derivative.
%                 Default: 20
%
%   'stepFunLen': One-sided length of the step function that is correlated
%                 with the time signal to obtain the temporal derivative
%                 signal, i.e. [-ones(1,stepFunLen) 0 ones(1,stepFunLen)].
%                 The step function represents a template for TMVs
%                 during depolarization.
%                 stepFunLen == 1 yields the central finite difference
%                 approximation of the first order derivative.
%                 Default: 1
%
%       'lambda': Parameter for spatial Laplacian smoothing.
%                 Spatial smoothing only affects the spatial derivative.
%                 Default: 1e2
%
%        'power': Power used for exponentiation of the derivative signal
%                 after normalization to a maximum value of 1.
%                 Increasing the power gives more weight to parts close to
%                 the maximum slope of the original signal.
%                 power == 1:      Classical correlation-based method.
%                 1 < power < inf: Correlation-based method with increased
%                                  emphasis on the intrinsic deflection.
%                 power == inf:    Purely deflection-based method
%                                  (no node pairs needed).
%                 Default: 1
%
%         'mesh': Mesh struct in the format of the vtkToolbox.
%                 Mandatory fields:
%                 mesh.points: [numPoints x 3] coordinates list,
%                 mesh.cells:  [numCells x 3] connectivity list.
%                 Can be omitted for power == inf AND derivative == 't'.
%
%        'pairs': Node pairs [numPairs x 2] used for cross-correlation.
%
% 'nodePairDist': Distance in number of edges for defining node pairs,
%                 if not directly provided via 'pairs'.
%                 Default: 2
%
%   'regression': Type of regression used to fit ATs to inter-node delays:
%                 'lad': least absolute deviations,
%                 'ls':  least squares.
%                 Default: 'lad'
%
%          'tol': Tolerance for iterative solver (linprog or lsqr).
%                 Default: 1e-8
%
%        'maxit': Maximum number of iterations.
%                 Default: 1000
%
% OUTPUT:
%             at: ATs in samples (1-based indexing) [numNodes x 1]
% 
% Copyright 2020 Steffen Schuler
% Institute of Biomedical Engineering
% Karlsruhe Institute of Technology
% www.ibt.kit.edu

%% Parse inputs

p = inputParser;
addParameter(p, 'upsampling', 10);
addParameter(p, 'derivative', 't');
addParameter(p, 'sigma', 20);
addParameter(p, 'stepFunLen', 1);
addParameter(p, 'lambda', 1e2);
addParameter(p, 'power', 1);
addParameter(p, 'mesh', []);
addParameter(p, 'pairs', []);
addParameter(p, 'nodePairDist', 2);
addParameter(p, 'regression', 'lad');
addParameter(p, 'tol', 1e-8);
addParameter(p, 'maxit', 1000);
parse(p, varargin{:});
p = p.Results;

if isempty(p.mesh)
    if contains(p.derivative, 's')
        error('Parameter ''mesh'' required to compute spatial derivative.');
    end
    if p.power ~= inf
        error('Parameter ''mesh'' required to compute node pairs.');
    end
end

%% Upsample and filter signal

sig = interp1(1:size(X,2), X', 1:1/p.upsampling:size(X,2))';
filtSig = gaussFiltfilt(sig, p.upsampling*p.sigma);

%% Define derivative signal

switch p.derivative
    case 's'
        derivSig = spatialDeriv(sig, p.mesh, p.lambda);
    case 't'
        derivSig = temporalDeriv(filtSig, p.stepFunLen);
    case 'st'
        derivSig = spatialDeriv(sig, p.mesh, p.lambda) .* temporalDeriv(filtSig, p.stepFunLen);
    otherwise
        error('Unknown derivative ''%s''.', p.derivative);
end

%% Actual activation time estimation

if p.power == inf 
    %% Deflection-based
    
    [~,at] = max(derivSig,[],2);
    at = (at-1)/p.upsampling+1;
else
    %% Correlation-based
    
    % Normalize and exponentiate derivSig
    derivSig = derivSig./repmat(max(derivSig,[],2),1,size(derivSig,2));
    derivSig = max(derivSig,0).^p.power;
    % for i = 1:size(derivSig,1)
    %     plot(derivSig(i,:))
    %     waitforbuttonpress
    % end
    
    % Compute node pairs
    if isempty(p.pairs)
        pairs = nodePairs(p.mesh, p.nodePairDist);
    else
        pairs = p.pairs;
    end
    
    % Compute delays
    delays = NaN(size(pairs,1),1);
    parfor i = 1:size(pairs,1)
        [xc,lag] = xcorr(derivSig(pairs(i,1),:), derivSig(pairs(i,2),:)); %#ok
        [~,ind] = max(xc);
        delays(i) = lag(ind);
    end
    delays = [delays; 0];
    
    % Compute difference matrix M
    numPoints = size(p.mesh.points,1);
    i = [repmat((1:size(pairs,1))',2,1); repmat(size(pairs,1)+1,numPoints,1)];
    j = [pairs(:,1); pairs(:,2); (1:numPoints)'];
    v = [ones(size(pairs,1),1); -ones(size(pairs,1),1); ones(numPoints,1)];
    M = sparse(i, j, v, size(pairs,1)+1, numPoints);

    % Compute activation times from delays using regression
    switch p.regression
        case 'lad'
            [at,flag,output] = lad(M, delays, p.tol, p.maxit);
            if ~flag
                warning('linprog failed with the following message:\n%s', output.message);
            end
        case 'ls'
            [at,flag] = lsqr(M, delays, p.tol, p.maxit);
            if flag
                warning('lsqr failed with flag %i.\n', flag);
            end
        otherwise
            error('Unknown regression ''%s''.', p.regression);
    end

    % Determine constant offset of activation times from mean of aligned signals
    derivSigAligned = delayseq(derivSig', -at);
    [~,offset] = max(mean(derivSigAligned,2));
    at = (at+offset-1)/p.upsampling+1;
end

end

%% Function definitions

function derivSig = temporalDeriv(sig, len)

switch len
    case 1
        s = [sig(:,1) sig sig(:,end)];
        derivSig = s(:,3:end)-s(:,1:end-2);
    otherwise
        template = [-ones(1,len) 0 ones(1,len)]; % must be mean-free!
        s = [repmat(sig(:,1),1,len) sig repmat(sig(:,end),1,len)];
        n = size(s,2);
        derivSig = NaN(size(sig));
        for i = 1:size(s,1)
            r = xcorr(s(i,:), template);
            derivSig(i,:) = r(n:end-2*len);
        end
end

end

function derivSig = spatialDeriv(sig, mesh, lambda)

if lambda > 0
    L = Laplacian(mesh);
    I = speye(size(L));
    sig = (I+lambda*(L'*L))\sig;
end

[~,Gx,Gy,Gz] = Gradient(mesh);
derivSig = sqrt((Gx*sig).^2+(Gy*sig).^2+(Gz*sig).^2);

end

function pairs = nodePairs(mesh, nodePairDist)

TR = vtkToTriangulation(mesh);
edg = TR.edges;
dist = distances(graph(edg(:,1), edg(:,2)));

numPoints = size(mesh.points,1);
pairs = NaN(20*numPoints,2);
numPairs = 0;
for i = 1:numPoints
    neighs = find(dist(:,i) == nodePairDist);
    numNeighs = numel(neighs);
    pairs(numPairs+1:numPairs+numNeighs,:) = [repmat(i,numNeighs,1) neighs];
    numPairs = numPairs + numNeighs;
end
pairs(isnan(pairs(:,1)),:) = [];
pairs = unique(sort(pairs,2), 'rows');

end

function [x,flag,output] = lad(A, b, tol, maxit)

c = [zeros(size(A,2),1); ones(size(b))];
F = [A -speye(size(b,1)); -A -speye(size(b,1))];
g = [b; -b];
options = optimoptions('linprog', 'Algorithm','interior-point', 'Display','off', 'OptimalityTolerance',tol, 'MaxIterations',maxit);
[z,~,flag,output] = linprog(c, F, g, [], [], [], [], options);
x = z(1:size(A,2));

end

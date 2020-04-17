function [G,Gx,Gy,Gz] = Gradient(vtk, atPoints)
% Relies on the gptoolbox by Alec Jacobson

if nargin < 2
    atPoints = true;
end

P = double(vtk.points);
C = double(vtk.cells);
G = grad(P,C);

if atPoints
    % compute area/volume weighted mean of cell values adjacent to a point
    i = C(:);
    j = repmat((1:size(C,1))', size(C,2), 1);
    if size(C,2)==4
        v = volume(P,C);
    elseif size(C,2)==3
        v = doublearea(P,C)/2;
    else
        error('Invalid number of columns in cell list.');
    end
    v = repmat(v, size(C,2), 1);
    W = sparse(i, j, v, size(P,1), size(C,1));
    W = spdiags(1./sum(W,2), 0, size(W,1), size(W,1)) * W;
    G = reshape(W*reshape(G, size(C,1), []), [], size(P,1));
end

if nargout > 1
    n = size(G,1)/3;
    Gx = G(1:n,:);
    Gy = G(1+n:2*n,:);
    Gz = G(1+2*n:end,:);
end

end
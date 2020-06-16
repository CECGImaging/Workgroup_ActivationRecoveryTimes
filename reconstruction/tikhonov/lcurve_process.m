function [optLambda, etaSmoothed, curvature, dcurvature] = lcurve_process(rho, eta, lambdas, smoothingParam, ax1, ax2, ax3)

if nargin < 4 || isempty(smoothingParam)
    smoothingParam = 1;
end

logRho = log10(rho);
logEta = log10(eta);
ft = fittype('smoothingspline');
opts = fitoptions('Method', 'SmoothingSpline');
opts.SmoothingParam = smoothingParam;

lcurvefit = fit(logRho, logEta, ft, opts);
logEtaSmoothed = feval(lcurvefit, logRho);
etaSmoothed = 10.^logEtaSmoothed;
[d1,d2] = differentiate(lcurvefit, logRho);
curvature = d2./(1+d1.^2).^(3/2);

curvaturefit = fit(logRho, curvature, ft, opts);
dcurvature = differentiate(curvaturefit, logRho);

[~,optInd] = findpeaks(dcurvature, 'MinPeakHeight', 0.1*max(dcurvature), 'NPeaks',1);
optLambda = lambdas(optInd);

if nargin > 4 && ~isempty(ax1)
    k = round(linspace(1, numel(lambdas), 10));
    loglog(ax1, rho, eta, 'b', rho, etaSmoothed, 'r', rho(k), etaSmoothed(k), 'k.', 'MarkerSize',10);
    text(ax1, rho(k), etaSmoothed(k), sprintfc(' %.2e', lambdas(k)));
    xlabel(ax1, 'Residual norm'); ylabel(ax1, 'Solution norm'); title(ax1, 'L-curve');
end
if nargin > 5 && ~isempty(ax2)
    [~,k] = findpeaks(curvature);
    semilogx(ax2, rho, curvature, 'r', rho(k), curvature(k), 'k.', 'MarkerSize',10);
    text(ax2, rho(k), curvature(k), sprintfc(' %.2e', lambdas(k)));
    xlabel(ax2, 'Residual norm'); ylabel(ax2, 'Curvature'); title(ax2, 'Curvature');
end
if nargin > 6 && ~isempty(ax3)
    semilogx(ax3, rho, dcurvature, 'r', rho(optInd), dcurvature(optInd), 'k.', 'MarkerSize',10);
    text(ax3, rho(optInd), dcurvature(optInd), sprintfc(' %.2e', optLambda));
    xlabel(ax3, 'Residual norm'); ylabel(ax3, 'Derivative of curvature'); title(ax3, 'Derivative of curvature');
end

end
addpath(genpath('..'));

%% Load input

% heart = vtkRead('geometry/heart_res1.ply');
load('geometry/heart_res1.mat');

X_true = load('signals/lvLateral/tmv_res1.mat');
X_true = double(X_true.tmv);
rng(1);
X = addwhitenoise(X_true, 0);

%% AT estimation

upsampling = 10;
sigma = 10;
lambda = 1e2;

at_true    = at_unified(X_true, 'upsampling',upsampling, 'sigma',1, 'power',inf);
at_defl_t  = at_unified(X, 'upsampling',upsampling, 'sigma',sigma, 'power',inf);
at_defl_st = at_unified(X, 'upsampling',upsampling, 'sigma',sigma, 'power',inf, 'derivative','st', 'mesh',heart, 'lambda',lambda);
at_stepFun = at_unified(X, 'upsampling',upsampling, 'sigma',sigma, 'power',inf, 'stepFunLen',300);
at_corr_t  = at_unified(X, 'upsampling',upsampling, 'sigma',sigma, 'power',1, 'mesh',heart, 'nodePairDist',2, 'regression','lad');
at_corr_st = at_unified(X, 'upsampling',upsampling, 'sigma',sigma, 'power',1, 'derivative','st', 'mesh',heart, 'lambda',lambda, 'nodePairDist',2, 'regression','lad');

%% Visualization

limits = [min(at_true) max(at_true)];
numSteps = 20;
angles = [50 60];

im_true    = atImage(heart, at_true,    limits, numSteps, angles, 'Truth');
im_defl_t  = atImage(heart, at_defl_t,  limits, numSteps, angles, sprintf('Defl T    CC = %.3f', corr(at_true, at_defl_t)));
im_defl_st = atImage(heart, at_defl_st, limits, numSteps, angles, sprintf('Defl ST    CC = %.3f', corr(at_true, at_defl_st)));
im_stepFun = atImage(heart, at_stepFun, limits, numSteps, angles, sprintf('StepFun    CC = %.3f', corr(at_true, at_stepFun)));
im_corr_t  = atImage(heart, at_corr_t,  limits, numSteps, angles, sprintf('Corr T    CC = %.3f', corr(at_true, at_corr_t)));
im_corr_st = atImage(heart, at_corr_st, limits, numSteps, angles, sprintf('Corr ST    CC = %.3f', corr(at_true, at_corr_st)));

blank = repmat(uint8(255), size(im_true));
im = [im_true im_defl_t im_defl_st im_stepFun; blank im_corr_t im_corr_st blank];

imshow(im);
imwrite(im, 'example_whiteNoise.png');

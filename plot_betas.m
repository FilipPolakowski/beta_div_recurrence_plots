% Author: Elena Deckert
% Contributors: Dr. ir. Martijn Boussé
% Version: Version 1.0 - 2024-14-06

function [R_1, R_2, R_3, R_4, R_5] = plot_betas(x, m, t, varargin)
% Computes 5 beta divergence RPs for beta = 0,0.5,1.1.5,2 and plots these.
%
% Inputs:
%   x: Time series.
%   m: Embedding dimension.
%   t: Time delay.
%   threshold: (Optional) Boolean deciding if RP should be thresholded or not. Default false.
%   epsilon: (Optional) Threshold parameter. Default 0.1.
%
% Output:
%   R: Recurrence matrices for 5 different beta values (0,0.5,1.1.5,2).
%% input parser
p = inputParser;

addRequired(p, 'x', @isnumeric);
addRequired(p, 'm', @isnumeric);
addRequired(p, 't', @isnumeric);

addOptional(p, 'threshold', false, @islogical); % default false
addOptional(p, 'epsilon', 0.1, @isnumeric); % default 0.1

%% initialize
parse(p,x,m,t, varargin{:});

x = p.Results.x;
m = p.Results.m;
t = p.Results.t;
threshold = p.Results.threshold;
epsilon = p.Results.epsilon;

%% calc beta recurrence matrices
R_1 = computeRM_blockedkronecker_nonsymmetric(x, m, t, 50, 0);
R_2 = computeRM_blockedkronecker_nonsymmetric(x, m, t, 50, 0.5);
R_3 = computeRM_blockedkronecker_nonsymmetric(x, m, t, 50, 1);
R_4 = computeRM_blockedkronecker_nonsymmetric(x, m, t, 50, 1.5);
R_5 = computeRM_blockedkronecker_nonsymmetric(x, m, t, 50, 2);

%% plot figure
figure('Position', [100, 42, 1000, 600]);

% use thresholded R if requested
if threshold
    R_1 = R_1 > epsilon;
    R_2 = R_2 > epsilon;
    R_3 = R_3 > epsilon;
    R_4 = R_4 > epsilon;
    R_5 = R_5 > epsilon;
    colormap([0 0 0 ; 1 1 1]);
end

%plot
subplot(2,3,1);
imagesc(R_1);
xlabel('Time');
ylabel('Time');
axis xy
axis tight
axis square

subplot(2,3,2);
imagesc(R_2);
xlabel('Time');
ylabel('Time');
axis xy
axis tight
axis square

subplot(2,3,3);
imagesc(R_3);
xlabel('Time');
ylabel('Time');
axis xy
axis tight
axis square

subplot(2,3,4);
imagesc(R_4);
xlabel('Time');
ylabel('Time');
axis xy
axis tight
axis square

subplot(2,3,5);
imagesc(R_5);
xlabel('Time');
ylabel('Time');
axis xy
axis tight
axis square

subplot(2,3,6);
plot(x);
xlabel('Time');
ylabel('Amplitude');
axis xy
axis tight
axis square


colorbar('Location', 'east', 'Position', [0.92, 0.1, 0.02, 0.8]);

title(subplot(2, 3, 1), 'beta = 0');
title(subplot(2, 3, 2), 'beta = 0.5');
title(subplot(2, 3, 3), 'beta = 1');
title(subplot(2, 3, 4), 'beta = 1.5');
title(subplot(2, 3, 5), 'beta = 2');
title(subplot(2, 3, 6), 'time series');

end
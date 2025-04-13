% Author: Elena Deckert
% Contributors: Dr. ir. Martijn Boussé
% Version:Version 1.0 - 2024-14-02

function [plot] = rp_plot(R, varargin)
% Plots the recurrence matrix
%
% Inputs:
%   R: Recurrence matrix.
%   epsilon: (Optional) Threshold parameter. Default unthresholded.
%
% Output:
%   plot: Recurrence plot.

%% parser options
p = inputParser;

addRequired(p, 'R', @isnumeric);

addOptional(p, 'epsilon',[], @isnumeric); %default []

parse(p, R, varargin{:});

%% initialsie variables
R = p.Results.R;
epsilon = p.Results.epsilon;

%% plot the function
 % if we do not have a threshold plot the values with  colorbar
 if isempty(epsilon)
    figure;
    imagesc(R);
    colorbar
 % if we have a threshold create binary recurrence matrix and plot lack and white
 else
    plot = R > epsilon; % create binary recurrence matrix
    figure;
    imagesc(plot);
    colormap([0 0 0 ; 1 1 1]); % White for non-recurrent, Black for recurrent
 end
    
 xlabel('Time');
 ylabel('Time');

% Adjust axis properties
axis xy
axis tight

end
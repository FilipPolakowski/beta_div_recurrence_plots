% Author: Elena Deckert
% Contributors: Dr. ir. Martijn Boussé
% Version:Version 1.0 - 2024-14-02

function [R] = rp_create(x, m, t, varargin)
% Calls the computation and plotting function to generate the RP
%
% Inputs:
%   x: Time series.
%   m: Embedding dimension.
%   t: Time delay.
%   metric: (Optional) Calculate euclidean distance ('norm') or beta divergence ('betadiv'). Default 'norm'.
%   beta: (Optional) Parameter controlling the divergence calculation. Default is [].
%   Kronecker: (Optional) Boolean deciding if Kronecker computation or loop. Default true.
%   epsilon: (Optional) Threshold parameter. Default [].
%   plot: (Optional) Boolean deciding RP is going to be plotted. Default true.
%
% Output:
%   R: Recurrence matrix.

%% parser options
    p = inputParser;

    addRequired(p, 'x', @isnumeric);
    addRequired(p, 'm', @isnumeric);
    addRequired(p, 't', @isnumeric);

    addOptional(p, 'epsilon', [], @isnumeric); % default []
    addOptional(p, 'metric', 'norm', @(x) ischar(x) && (strcmpi(x, 'norm') || strcmpi(x, 'betadiv'))); % default norm
    addOptional(p, 'Kronecker', true, @islogical); % default true
    addOptional(p, 'beta',[], @isnumeric); % default []
    addOptional(p, 'plot', true, @islogical); % default true


    parse(p, x, m, t, varargin{:});

%% initialise variables
    x = p.Results.x;
    m = p.Results.m;
    t = p.Results.t;
    epsilon = p.Results.epsilon;
    metric = p.Results.metric;
    Kronecker = p.Results.Kronecker;
    beta = p.Results.beta;
    plot = p.Results.plot;


%% input testing
    input_test(x,m,t,'epsilon', epsilon);

%% call computation function
    % decide if kronecker computation or not 
    if Kronecker == false  
       R = rp_compute(x, m, t, 'metric', metric, 'beta', beta);
    else
       R = rp_kronecker_compute(x, m, t, 'metric', metric, 'beta', beta);
    end
    
%% call plot function
    if plot
        rp_plot(R, 'epsilon', epsilon);
    end
   

end
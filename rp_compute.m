%Author: Elena Deckert
% Contributors: Dr. ir. Martijn Boussé
% Version:Version 1.0 - 2024-14-02

function [R] = rp_compute(x, m, t, varargin)
% Conducts the recurrence matrix computations with a nested loop
%
% Inputs:
%   x: Time series.
%   m: Embedding dimension.
%   t: Time delay.
%   metric: (Optional) Calculate euclidean distance ('norm') or beta divergence ('betadiv'). Default 'norm'.
%   beta: (Optional) Parameter controlling the divergence calculation. Default is [].
%
% Output:
%   R: Recurrence matrix.

%% parser options
p = inputParser;

addRequired(p, 'x', @isnumeric);
addRequired(p, 'm', @isnumeric);
addRequired(p, 't', @isnumeric);

addOptional(p, 'metric', 'norm', @(x) ischar(x) && (strcmpi(x, 'norm') || strcmpi(x, 'betadiv')));
addOptional(p, 'beta',[], @isnumeric);

parse(p, x, m, t, varargin{:});

%% initialise variables

x = p.Results.x;
m = p.Results.m;
t = p.Results.t;
metric = p.Results.metric;
beta = p.Results.beta;

N = length(x);
n = N-(m-1)*t; 
R = zeros(n, n);

%% Recurrence Matrix computation

for i = 1:n
    for j = 1:n
        
    % Embedding vectors: takes x from i-th to the (i + (m-1)*t)-th element in steps of t
        vec1 = x(i:t:i + (m-1)*t);
        vec2 = x(j:t:j + (m-1)*t);            
        
    % Distance Calculation
        % Euclidean distance
        if strcmpi(metric, 'norm') || isempty(metric)
            dist = norm(vec1 - vec2);
        % Beta Divergences
        elseif strcmpi(metric, 'betadiv')
            dist = beta_calc(vec1, vec2, 'beta', beta);
        else
            error('Invalid metric. Supported values are ''norm'' or ''betadiv''.');
        end
        
        % Update recurrence matrix
        R(i, j) = dist;
        
    end
end
R = abs(R); % use absolute value because the direction of distance is irrelevant
R = rescale(R); % rescale values between 0 and 1

end

    
    
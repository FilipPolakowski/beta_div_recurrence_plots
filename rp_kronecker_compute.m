%Author: Elena Deckert
% Contributors: Dr. ir. Martijn Boussé
% Version:Version 1.0 - 2024-14-02

function [R] = rp_kronecker_compute(x, m, t, varargin)
% Conducts the recurrence matrix computations with the kronecker product
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
H = zeros(n, m); % H matrix used for kronecker calculations
 
%% Recurrence Matrix computation
% create H matrix: stacking the embedding vectors 
    for i= 1:n
        vec = x(i:t:i + (m-1)*t);
        H(i,:) = vec;
    end

% Kronecker product 
    vec1 = kron(ones(n,1), H);
    vec2 = kron(H, ones(n,1));

% Calculate distance
    % Euclidean distance
    if strcmpi(metric, 'norm') || isempty(metric)
        dist = vecnorm(vec1 - vec2, 2, 2); %row wise norm calc
    % Beta divergence
    elseif strcmpi(metric, 'betadiv')
        dist = beta_calc(vec1, vec2, 'beta', beta);
    else
        error('Invalid metric. Supported values are ''norm'' or ''betadiv''.');
    end

R = reshape(dist, n, n);
R = abs(R); % use absolute value because the direction of distance is irrelevant
R = rescale(R); %rescale values between 0 and 1
end

    
% Author: Elena Deckert
% Contributors: Dr. ir. Martijn Boussé
% Version:Version 1.0 - 2024-14-02

function dist = beta_calc(vec1, vec2, varargin)
% Calculates the beta divergence between two input vectors.
%
% Inputs:
%   vec1: First embedding vector.
%   vec2: Second embedding vector.
%   beta: (Optional) Parameter controlling the divergence calculation. Default value is beat=2.
%
% Output:
%   dist: The calculated beta divergence value between vec1 and vec2.
%% parser options
p = inputParser;

addRequired(p, 'vec1', @isnumeric);
addRequired(p, 'vec2', @isnumeric);

addOptional(p, 'beta',2, @isnumeric); % default value 2 (euclidean)

parse(p, vec1, vec2, varargin{:});

%% initialise variables
vec1 = p.Results.vec1;
vec2 = p.Results.vec2;  
beta = p.Results.beta;

%% calculate beta divergence
     one = ones(size(vec1)); %create a vector of ones to use for vector calculations

     if beta == 0 %IS
        beta_vector = (vec1 ./ vec2)-log(vec1 ./ vec2) - one;

    elseif beta == 1 %KL
        beta_vector = vec1 .* (log(vec1)-log(vec2)) + vec1 - vec2;

     else 
        beta = beta.*one; %create a vector of the value of beta to use for vector calculations
        beta_vector = (vec1.^beta + (beta-one) .* vec2.^beta - beta .* vec1 .* vec2.^(beta-one)) ./ (beta .* (beta-one));

    end
    dist = sum(beta_vector,2); % sum over all single divergence values to get vector divergence
end
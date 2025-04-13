% Author: Elena Deckert
% Contributors:Dr. ir. Martijn Boussé
% Version: Version 1.0 - 2024-14-06

function dist = beta_calc_e(vec1, vec2, beta)
% Calculates the beta divergence between two input vectors.
% Faster than beta_calc because there is no optional parameters
%
% Inputs:
%   vec1: First embedding vector.
%   vec2: Second embedding vector.
%   beta: Controlling the divergence calculation.
%
% Output:
%   dist: The calculated beta divergence value between vec1 and vec2.

%% calculate beta divergence
     one = ones(size(vec1)); %create a vector of ones to use for vector calculations

     if beta == 0 % IS
        temp = vec1 ./ vec2;
        beta_vector = (temp)-log(temp) - one;

    elseif beta == 1 % KL
        beta_vector = vec1 .* (log(vec1)-log(vec2)) + vec1 - vec2;

     else 
        beta = beta.*one; %create a vector of the value of beta to use for vector calculations
        beta_vector = (vec1.^beta + (beta-one) .* vec2.^beta - beta .* vec1 .* vec2.^(beta-one)) ./ (beta .* (beta-one));

    end
    dist = sum(beta_vector,2); % sum over all single divergence values to get vector divergence
end
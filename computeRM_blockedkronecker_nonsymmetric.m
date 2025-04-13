
% Author: Svyatoslav Faymonville, Dr. Ir. Martijn Boussé, Elena Deckert
% Contributors:
% Version: version 1 (14.03.2024, Svety) kronecker but blocked into 4 blocks.
%          version 1.0.1 blocked arbitrarely
%          version 1.0.2 (02.04.2024) added input check
%          version 2 (13.06.2024) using beta divergences

function RM = computeRM_blockedkronecker_nonsymmetric(x,M,T,mS,beta)
% Computes the recurrence matrix with beta-divergences using blocking for faster computations.
%
% Inputs:
%   x: Time Series.
%   M: Embedding dimension.
%   T: Time delay.
%   mS: the number of segments the time sries is divided into (recommended mS > 10).
%
% Output:
%   RM: Recurrence Matrix.
%%
[n,m] = size(x);


%check what kind of vector: row or column
if n>m
    x = x';
    [n,m] = size(x);
end

% error message for wrong input dimensions
if n>1
    msg = 'wrong input dimensions';
    error(msg)
end

MT = (M - 1) * T;
X = zeros(m-MT,M);
for i = 1:m-MT
    X(i,:) = x( 1 , i : T : i + MT );
end

segments = cell(mS,1,1);
for cS = 1:mS
    start_point = 1 + floor((cS-1)/mS * (m-MT)) ;
    end_point = floor(cS/mS * (m-MT)) ;
    segments{cS} = X(start_point:end_point , : ) ;
end

blocks = cell(mS, mS) ;

% calculates all blocks
for i = 1:mS
    for j = 1:mS
        [Li,~] = size(segments{i});
        [Lj,~] = size(segments{j});

        Vi = ones(Li,1);
        Vj = ones(Lj,1);
        
        vec1 = kron(Vi,segments{j});
        vec2 = kron(segments{i},Vj);

        dist = beta_calc_e(vec1, vec2, beta); % added beta as distance metrix
        R = reshape(dist,Lj,Li);
        R = abs(R);

        blocks{j,i} = R;
    end
end

% puts all blocks together into the RM
for i = 1:mS
    blocks{i} = cat(2,blocks{i,:});
end
RM = cat(1,blocks{:,1});
RM = rescale(RM);

end
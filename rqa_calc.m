% Author: Elena Deckert
% Contributors: Dr. ir. Martijn Boussé
% Version:Version 1.0 - 2024-14-04
function [rr, det, l, lam, tt, hor, ah] = rqa_calc(R, varargin)
% Computes the RQA measures from the recurrence matrix
%
% Inputs:
%   R: Recurrence magtrix.
%   lmin: (Optional) Minimum line length. Default is 5.
%   epsilon: (Optional) Threshold parameter. Default 0.1.
%   rr, det, l, lam, tt: (Optional) Boolean deciding if RQA should be calculated or not. Default we calc all RQA measures.
%
% Output:
%   rr, det, l, lam, tt: values of RQA measures.
%% parser options
p = inputParser;

addRequired(p, 'R', @isnumeric);

addOptional(p, 'epsilon',0.1, @isnumeric); % default set to 0.1
addOptional(p, 'lmin', 5, @isnumeric); % default set to 5
addOptional(p, 'rr', true, @islogical); % default true
addOptional(p, 'det', true, @islogical);% default true
addOptional(p, 'l', true, @islogical);% default true
addOptional(p, 'lam', true, @islogical);% default true
addOptional(p, 'tt', true, @islogical);% default true
addOptional(p, 'hor', true, @islogical);% default true
addOptional(p, 'ah', true, @islogical);% default true

parse(p, R, varargin{:});

%% initialsie variables
R = p.Results.R;
rr = p.Results.rr;
det = p.Results.det;
lam = p.Results.lam;
l = p.Results.l;
tt = p.Results.tt;
hor = p.Results.hor;
ah = p.Results.ah;
epsilon = p.Results.epsilon;
lmin = p.Results.lmin;

R_thresholded = R < epsilon; % create thresholded RP
R_thresholded = rot90(R_thresholded); % rotate so it represents graph and is easier to understand calculations
[rows, cols] = size(R_thresholded);

%% calc rr
%recurrence rate is the percentage of recurrence points, so ones in the thresholded matrix
if rr
    rr = sum(R_thresholded(:)) / numel(R_thresholded);
end

%% calc det
%determinism is the proportion of recurrent points vs recurrent points in diagonals
if det 
    % Count the number of diagonal lines
    num_diagonal_lines = 0;
    diagonal_lines = zeros(1, rows*cols);
    
    % Traverse diagonally from top-left to bottom-right 
    for k = 1:(rows + cols - 1)

        % Skip the central diagonal
        if k == ceil((rows + cols) / 2)
            continue;
        end
        
        count = 0;
        % Determine starting row and column for this diagonal
        if k <= rows
            start_row = k;
            start_col = 1;
        else
            start_row = rows;
            start_col = k - rows + 1;
        end
        % Count consecutive 1 on this diagonal
        i = start_row;
        j = start_col;
        while i >= 1 && j <= cols
            if R_thresholded(i, j) == 1
                count = count + 1;
            else
                if count >= lmin % add count if its longer than lmin
                    num_diagonal_lines = num_diagonal_lines + 1;
                    diagonal_lines(num_diagonal_lines) = count;
                end
                count = 0;
            end
            i = i - 1;
            j = j + 1;
        end
        if count >= lmin
            num_diagonal_lines = num_diagonal_lines + 1;
            diagonal_lines(num_diagonal_lines) = count;
        end
    end
    
    valid_lines = diagonal_lines(diagonal_lines > 0);
    det = sum(valid_lines(:)) / (sum(R_thresholded(:)) - sum(diag(R_thresholded))); 
end

%% calc l
% is the average length of diagonal lines
if l 
    num_diagonal_lines = 0;
    diagonal_lines = zeros(1, rows*cols);
    
    % Traverse diagonally from top-left to bottom-right
    for k = 1:(rows + cols - 1)

        % Skip the central diagonal
        if k == ceil((rows + cols) / 2)
            continue;
        end

        count = 0;
        % Determine starting row and column for this diagonal
        if k <= rows
            start_row = k;
            start_col = 1;
        else
            start_row = rows;
            start_col = k - rows + 1;
        end
        % Count consecutive 1 on this diagonal
        i = start_row;
        j = start_col;
        while i >= 1 && j <= cols
            if R_thresholded(i, j) == 1
                count = count + 1;
            else
                if count >= lmin % add count if its longer than lmin
                    num_diagonal_lines = num_diagonal_lines + 1;
                    diagonal_lines(num_diagonal_lines) = count;
                end
                count = 0;
            end
            i = i - 1;
            j = j + 1;
        end
        if count >= lmin
            num_diagonal_lines = num_diagonal_lines + 1;
            diagonal_lines(num_diagonal_lines) = count;
        end
    end
    
    valid_lines = diagonal_lines(diagonal_lines > 0);
    if isempty(valid_lines)
        l = 0;
    else
        l = mean(valid_lines);
    end
end


%% calc lam
%laminarity is the proportion of points in vertical recurrent lines vs. recurrent points
if lam 
    num_vertical_lines = 0;
    vertical_lines = zeros(1, rows*cols);
    
    % traverse
    for j = 1:cols
        count = 0;
        for i = 1:rows
            if R_thresholded(i, j) == 1
                count = count + 1;
            else
                if count >= lmin % add count if its longer than lmin
                    num_vertical_lines = num_vertical_lines + 1;
                    vertical_lines(num_vertical_lines) = count;
                end
                count = 0;
            end
        end
        if count >= lmin
            num_vertical_lines = num_vertical_lines + 1;
            vertical_lines(num_vertical_lines) = count;
        end
    end
    
    valid_lines = vertical_lines(vertical_lines > 0);
    lam = sum(valid_lines(:)) / sum(R_thresholded(:));
end

%% calc tt
% average length of vertical lines
if tt 
    num_vertical_lines = 0;
    vertical_lines = zeros(1, rows*cols);
    
    % traverse
    for j = 1:cols
        count = 0;
        for i = 1:rows
            if R_thresholded(i, j) == 1
                count = count + 1;
            else
                if count >= lmin % add count if its longer than lmin
                    num_vertical_lines = num_vertical_lines + 1;
                    vertical_lines(num_vertical_lines) = count;
                end
                count = 0;
            end
        end
        if count >= lmin
            num_vertical_lines = num_vertical_lines + 1;
            vertical_lines(num_vertical_lines) = count;
        end
    end
    
    valid_lines = vertical_lines(vertical_lines > 0);
    
    if isempty(valid_lines)
        tt = 0;
    else
        tt = mean(valid_lines);
    end

end

%% calc hor
% is the proportion of points in horizontal recurrent lines vs. recurrent points
if hor 
    R_hor = rot90(R_thresholded);
    num_vertical_lines = 0;
    vertical_lines = zeros(1, rows*cols);
    
    % traverse
    for j = 1:cols
        count = 0;
        for i = 1:rows
            if R_hor(i, j) == 1
                count = count + 1;
            else
                if count >= lmin % add count if its longer than lmin
                    num_vertical_lines = num_vertical_lines + 1;
                    vertical_lines(num_vertical_lines) = count;
                end
                count = 0;
            end
        end
        if count >= lmin
            num_vertical_lines = num_vertical_lines + 1;
            vertical_lines(num_vertical_lines) = count;
        end
    end
    
    valid_lines = vertical_lines(vertical_lines > 0);
    hor = sum(valid_lines(:)) / sum(R_thresholded(:));
end

%% calc ah
% average length of horizontal lines
if ah
    R_ah = rot90(R_thresholded);
    num_vertical_lines = 0;
    vertical_lines = zeros(1, rows*cols);
    
    % traverse
    for j = 1:cols
        count = 0;
        for i = 1:rows
            if R_ah(i, j) == 1
                count = count + 1;
            else
                if count >= lmin % add count if its longer than lmin
                    num_vertical_lines = num_vertical_lines + 1;
                    vertical_lines(num_vertical_lines) = count;
                end
                count = 0;
            end
        end
        if count >= lmin
            num_vertical_lines = num_vertical_lines + 1;
            vertical_lines(num_vertical_lines) = count;
        end
    end
    
    valid_lines = vertical_lines(vertical_lines > 0);
    
    if isempty(valid_lines)
        ah = 0;
    else
        ah = mean(valid_lines);
    end

end

end
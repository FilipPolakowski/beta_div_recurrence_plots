% Author: Elena Deckert
% Contributors: Dr. ir. Martijn Boussé
% Version:Version 1.0 - 2024-14-06
function [rr, det, l, lam, tt, hor, ah, lmax, vmax, hmax, t2_v, t2_h, entr] = rqa_calc_e(R, lmin, epsilon)
% Computes the RQA measures from the recurrence matrix
% faster and more RQAs then rqa_calc
%
% Inputs:
%   R: Recurrence magtrix.
%   lmin: (Optional) Minimum line length. Default is 5.
%   epsilon: Threshold parameter. 
%
% Output:
%   rr, det, l, lam, tt, hor, ah, lmax, vmax, hmax, t2_v, t2_h, entr: values of RQA measures.

R_thresholded = R < epsilon; % create thresholded RP
R_thresholded = rot90(R_thresholded); % rotate so it represents graph and is easier to understand calculations
[rows, cols] = size(R_thresholded);

%% calc rr
%recurrence rate is the percentage of recurrence points, so ones in the thresholded matrix
    rr = sum(R_thresholded(:)) / numel(R_thresholded);

%% calc det, lmax, entr
%determinism is the proportion of recurrent points vs recurrent points in diagonals
%l is average diagonal line length
%lmax is the max diagonal line length
%entropy is the probability distribution of diagonal line lengths

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
            if R_thresholded(i, j) 
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
    
    %det
        det = sum(valid_lines(:)) / (sum(R_thresholded(:)) - sum(diag(R_thresholded))); 
    %l
        if isempty(valid_lines)
            l = 0;
        else
            l = mean(valid_lines);
        end
    
    %lmax
        lmax = max(diagonal_lines);
        
    %entropy
        if isempty(valid_lines)
            entr = 0;
        else 
            unique_lengths = unique(valid_lines);
            counts = histcounts(valid_lines, unique_lengths);
            p_l = counts / length(valid_lines);
            entr = -sum(p_l .* log(p_l));
        end

%% calc lam, tt, hor, ah, vmax, hmax, t2_v, t2_h
%laminarity is the proportion of points in vertical recurrent lines vs. recurrent points
%tt is the average vertical line length 
%vmax is the max vertical line length 
%hor is the proportion of points in horizontal recurrent lines vs. recurrent points
%ah is the average horizontal line length 
%hmax is the max horozontal line length 
%t2_v average length of white vertical lines
%t2_h average length of white horizontal lines

    num_vertical_lines = 0;
    vertical_lines = zeros(1, rows*cols);
    
    num_horizontal_lines = 0;
    horizontal_lines = zeros(1, rows*cols);
    
    num_vertical_lines_w = 0;
    vertical_lines_w = zeros(1, rows*cols);
    
    num_horizontal_lines_w = 0;
    horizontal_lines_w = zeros(1, rows*cols);

    % traverse
    for j = 1:cols
        count_v = 0;
        count_h = 0;
        count_white_v = 0;
        count_white_h = 0;
        
        for i = 1:rows
            %vertical
            if R_thresholded(i, j) 
                count_v = count_v + 1;
            else
                if count_v >= lmin % add count if its longer than lmin
                    num_vertical_lines = num_vertical_lines + 1;
                    vertical_lines(num_vertical_lines) = count_v;
                end
                count_v = 0;
            end

            %vertical white
            if ~R_thresholded(i, j) 
                count_white_v = count_white_v + 1;
            else
                if count_white_v >= lmin % add count if its longer than lmin
                    num_vertical_lines_w = num_vertical_lines_w + 1;
                    vertical_lines_w(num_vertical_lines_w) = count_white_v;
                end
                count_white_v = 0;
            end
            
            %horizontal
            if R_thresholded(j, i) 
                count_h = count_h + 1;
            else
                if count_h >= lmin % add count if its longer than lmin
                    num_horizontal_lines = num_horizontal_lines + 1;
                    horizontal_lines(num_horizontal_lines) = count_h;
                end
                count_h = 0;
            end

            %horizontal white
            if ~R_thresholded(j, i) 
                count_white_h = count_white_h + 1;
            else
                if count_white_h >= lmin % add count if its longer than lmin
                    num_horizontal_lines_w = num_horizontal_lines_w + 1;
                    horizontal_lines_w(num_horizontal_lines_w) = count_white_h;
                end
                count_white_h = 0;
            end

        end

        if count_v >= lmin
            num_vertical_lines = num_vertical_lines + 1;
            vertical_lines(num_vertical_lines) = count_v;
        end

        if count_h >= lmin
            num_horizontal_lines = num_horizontal_lines + 1;
            horizontal_lines(num_horizontal_lines) = count_h;
        end

        if count_white_v >= lmin
            num_vertical_lines_w = num_vertical_lines_w + 1;
            vertical_lines_w(num_vertical_lines_w) = count_white_v;
        end

        if count_white_h >= lmin
            num_horizontal_lines_w = num_horizontal_lines_w + 1;
            horizontal_lines_w(num_horizontal_lines_w) = count_white_h;
        end
    end
    
    valid_lines_v = vertical_lines(vertical_lines > 0);
    lam = sum(valid_lines_v(:)) / sum(R_thresholded(:));
     if isempty(valid_lines_v)
        tt = 0;
    else
        tt = mean(valid_lines_v);
    end

    valid_lines_h = horizontal_lines(horizontal_lines > 0);
    hor = sum(valid_lines_h(:)) / sum(R_thresholded(:));
    if isempty(valid_lines_h)
        ah = 0;
    else
        ah = mean(valid_lines_h);
    end

    vmax = max(vertical_lines);
    hmax = max(horizontal_lines);

    valid_lines_v_w = vertical_lines_w(vertical_lines_w > 0);
     if isempty(valid_lines_v_w)
        t2_v = 0;
    else
        t2_v = mean(valid_lines_v_w);
    end

    valid_lines_h_w = horizontal_lines_w(horizontal_lines_w > 0);
    if isempty(valid_lines_h_w)
        t2_h = 0;
    else
        t2_h = mean(valid_lines_h_w);
    end

end
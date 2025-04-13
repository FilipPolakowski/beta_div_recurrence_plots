%Author: Elena Deckert
% Contributors: Dr. ir. Martijn Boussé
% Version:Version 1.0 - 2024-14-06

%% set parameters
m = 10;                   % Embedding dimension
t = 1;                    % Time delay
epsilon = 0.1;            % Threshold (optional)
beta = 1;                 % Beta value (optional)
signal = 0;              
% choose the signal:
% 1: sin
% 2: sum of sins
% 3: damped sin
% 4: sin with trend
% 5: simple step
%6: simple 2 pyramid
%7: two peaks
%8: two main peaks and little ones
%9: preprocessed audio

%% signals
switch signal
    case 0
        y = linspace(0, 5, 100);  % 100 points in the interval [0, 5]
        x = sin(2 * pi * y); %sin
    case 1
        y = linspace(0, 5, 100);  % 100 points in the interval [0, 5]
        x = sin(2 * pi * y) + 0.5 .* sin(pi * y) + 0.75 .* cos(10 * pi * y); %sum of sin
    case 2
        y = linspace(0, 5, 100);  % 100 points in the interval [0, 5]
        x = exp(-0.7*y) .* sin(2 * pi * y); %damped sin
    case 3
        y = linspace(0, 5, 100);  % 100 points in the interval [0, 5]
        x = 0.1 .* y + sin(2 * pi * y); %sin with trend
    case 4 
        y = linspace(0, 5, 100);  % 100 points in the interval [0, 5]
        x = randn(100,1); %random
    case 5 %simple step
        m = 3;
        %x = [6,6,6,6,6,2,2,2,2,2];
        x = [3,3,3,3,3,2,2,2,2,2];
    case 6
        m = 3;
        x = [1, 2, 3, 2, 1, 2, 3, 2, 1];
    case 7 
        m = 3;
        %x = [1,1,1,1,1,9,14,9,1,1,1,1,1,1,1];
        x = [1,1,1,1,9,1,1,1,1,4,1,1,1,1,1,1,1,1,1,1];
    case 8
        m = 3;
        x = [1,1,1,1,3,9,1,1,1,1,1,1,1,2,4,1,1,1,1,3,2,12,1,1,1,3,1,1,1,1,2,1,1,1,2,1,1,1,4,1,1,1,2,1,1,1,1,1,2,1];
        %x = x+50;
        %x = x.*50;
    case 9
        [audio, Fs] = audioread('StutterTalk_2_19.wav');
                
                newFs = 2000; 
                audio = resample(audio, newFs, Fs);
            
                x = audio.^2;
end

%% create recurrence plot

%R = rp_create(audio_signal, m, t, 'metric', 'norm', 'Kronecker', true, 'plot', true);
%R_beta = rp_create(x, 25, 1, 'metric', 'betadiv', 'Kronecker', true, 'beta', 1, 'plot', true, 'epsilon', 0.1);

%[rr, det, l, lam, tt] = rqa_calc(R);
%[rrb, detb, lb, lamb, ttb] = rqa_calc(R_beta);

plot_betas(x,m,t,'threshold', true, 'epsilon', epsilon);


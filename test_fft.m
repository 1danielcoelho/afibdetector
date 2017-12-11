clc
close all
clear all

Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
X = S + 2*randn(size(t));

figure;
plot(1000*t(1:50),X(1:50))
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('t (milliseconds)')
ylabel('X(t)')

% Using this formula for 'n' speeds up FFT calculations, since
% powers of two are preferred
n = 2^nextpow2(L);
Y = fft(X, n);

% Frequency domain for the single-sided spectrum will go from
% 0 to Fs/2, with half the points as n plus 1
f = Fs*(0:(n/2))/n;

% Normalize by dividing by Fs according to Parseval's theorem
% https://math.stackexchange.com/questions/636847/understanding-fourier-transform-example-in-matlab
P = abs(Y/Fs);

% Multiply by two since we took half the power away when we
% only took half the points. The DC power should stay the same though
% https://www.12000.org/my_notes/on_scaling_factor_for_ftt_in_matlab/
P1 = 2 * P(1:n/2+1);
P1(1) = P1(1) / 2;

figure;
plot(f, P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
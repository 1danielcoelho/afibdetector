clc
close all
clear all

rng default

fs = 1000;
t = 0:1/fs:5-1/fs;
x = cos(2*pi*100*t)+randn(size(t));

% According to https://www.mathworks.com/help/signal/ref/pwelch.html
% By default, x is divided into the longest possible segments
% to obtain as close to but not exceed 8 segments with 50% overlap

% Each segment is windowed with a Hamming window

% The modified periodograms are averaged to obtain the PSD estimate

% The default nfft is the greater of 256 or the next power of 2
% greater than the length of the segments. In our case this will be
% 1000/8 -> 125 -> 128 + DC -> 129 points
[pxx,f] = pwelch(x, [], [], [], fs);

figure;
plot(f,10*log10(pxx))

xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
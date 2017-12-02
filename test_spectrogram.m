% Clean workspace
clear all
clc
close all

% Create cosine
Ts = 0.001;
t = 0:Ts:5-0.001;
freq = 100;
s = cos(2*pi*freq*t);

% Add some noise so we can see differences in the columns
s = s + 0.01 * randn(1, length(s));

figure;
windowSize = 500;  % Chunk size
overlap = [];  % Samples to overlap, default is 50%
nfft = 128;  % How many FFT points desired per chunk
samplingFreq = 1 / Ts;
spectrogram(s, windowSize, overlap, nfft, samplingFreq, 'yaxis');


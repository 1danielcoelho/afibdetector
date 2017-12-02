clear all
clc
close all

n = 0:999;
x1 = cos(pi/4*n)+0.5*sin(pi/2*n);
fs = 1000;
t = 0:0.001:1-0.001;
x2 = cos(2*pi*125*t)+0.5*sin(2*pi*250*t);
freq = [125 250];

figure,
spectrogram(x2,[],[],[],fs,'yaxis')
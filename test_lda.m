clear all
clc
close all

% Load some sample data
load cities

% Normalize by their variance
w = 1./var(ratings);

[wcoeff,score,latent,tsquared,explained] = pca(ratings,...
'VariableWeights',w);

% Principal components are the columns of wcoeff. How much these
% components explain the data is in 'explained'
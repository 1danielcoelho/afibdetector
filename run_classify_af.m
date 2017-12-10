%% Cleanup
clear all
clc
close all

%% Run tests
trainingSetSizes = [0.1, 0.3, 0.5, 0.7, 0.9];
pcaPrincipalComponents = [0, 1, 3, 5, 10, 30, 50, 129];

for ts=trainingSetSizes
    for pc=pcaPrincipalComponents        
        classify_af('linear', ts, 99, pc)
    end
end
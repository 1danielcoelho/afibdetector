%% Cleanup
clc
clear all
close all

%% Import data
disp('Importing data...');

% Adds the WFDB toolkit folder to path so that we can use rdsamp and rdann
fullpath = mfilename('fullpath');
k = strfind(fullpath, mfilename);
fullfolder = fullpath(1:k-1);
addpath(genpath(fullfolder));

% Move to outside the WFDB folder so that it can find our records
prev_folder = pwd;
cd(fullfolder)

% Names of all the records we're going to use
recordNames = {'04015', '04043', '04048', '04126', '04746', '04908', ...
               '04936', '05091', '05121', '05261', '06426', '06453', ...
               '06995', '07162', '07859', '07879', '07910', '08215', ...
               '08219', '08378', '08405', '08434', '08455'};

learningSetCount = 0;
           
% Loads all the records into a 'records' structure
for recIndex = 1:length(recordNames)
    recName = char(recordNames(recIndex));
    disp(recName);
    recPath = strcat('mit-bih\', recName);    
    
    % Actually read everything    
    [r.signalmV,r.Fs,r.tmSecs]=rdsamp(recPath);
    [r.annSamples,r.anntype,r.subtype,r.chan,r.num,r.comments] = rdann(recPath, 'atr');
    r.annVec = get_af_annotation_vector(length(r.signalmV), r.annSamples, r.comments);
    
    % Mark it as a learning set
    r.isLearningSet = recIndex < 4;
    
    % Count the number of learning sets we have so far
    if r.isLearningSet
        learningSetCount = learningSetCount + 1;
    end
    
    % Pack everything we read into our records structure
    records.(strcat('rec', recName)) = r;    
end

% Go back to where we were before, if it matters
cd(prev_folder);

%% Break ECG1 signal records into windows
disp('Separating ECG1 into windows and extracting their classes...');
windowSizeSeconds = 4;

recordNames = fieldnames(records);
for recordName=recordNames'    
    record = records.(recordName{1});
    
    samplingFreq = record.Fs;
    windowSize = windowSizeSeconds * samplingFreq;
    numWindows = floor(length(record.signalmV(:, 1)) / windowSize);

    classes = zeros(numWindows, 1);
    windows = zeros(numWindows, windowSize);

    % Break record into windows
    for i=1:numWindows
        rangeStart = 1 + (i-1)*windowSize;
        rangeEnd = i*windowSize;    

        windows(i, :) = record.signalmV(rangeStart:rangeEnd, 1);    
        sampleAnns = record.annVec(rangeStart:rangeEnd);

        % If most of the samples are marked as AFIB, mark the window as AFIB
        classes(i) = sum(sampleAnns) > windowSize/2;
    end    
    
    records.(recordName{1}).signalmVWindows = windows;
    records.(recordName{1}).signalmVWindowsClasses = classes;
end

%% Get Welch PSD estimator from the windows and store them back into records
disp('Getting Welch PSD estimator from ECG1 windows...');

recordNames = fieldnames(records);
for recordName=recordNames'     
    record = records.(recordName{1});
    
    [psds, f] = pwelch(record.signalmVWindows',[],[],[],250);
    psds = psds';
    
    % Keep track of how many frequency bands pwelch returned. This has
    % been determined to be the same for all datasets
    numberFrequencyBands = size(psds, 2);
    
    %Store PSDs with each row being a window, each column being a frequency
    %(that is, a variable)
    %psds+(psds==0) first sets to 1 elements that are zero
    records.(recordName{1}).PSDs = 10*log10(psds+(psds==0));    
    
    nans = sum(isnan(records.(recordName{1}).PSDs));
    infs = sum(isinf(records.(recordName{1}).PSDs));
    
    %min(min(records.(recordName{1}).PSDs))
    
    if nans
        disp(strcat(recordName{1}, ' has NaNs in its PSDs'));
    end
    
    if infs
        disp(strcat(recordName{1}, ' has Infs in its PSDs'));
    end
end

%% Extract the learning sets from the data to train our SVM
disp('Separating learning set...');

trainingPSDs = zeros(1, numberFrequencyBands);
trainingClasses = zeros(1, 1);

recordNames = fieldnames(records);
for recordName=recordNames'     
    record = records.(recordName{1});
    
    if record.isLearningSet
        % Stack all training PSDs on the same matrix
        trainingPSDs = [trainingPSDs; record.PSDs];
        trainingClasses = [trainingClasses; record.signalmVWindowsClasses];        
    end
end

trainingPSDs = trainingPSDs(2:end, :);
trainingClasses = trainingClasses(2:end, :);

disp('Training set dimensions (samples x frequencies): ');
size(trainingPSDs)

%% Perform PCA of psds
disp('Performing PCA from window PSDs...');

numPrincipalComponents = 10;

[coeff,score,latent,tsquared,explained] = pca(trainingPSDs);

% Get first 10 components
reducedTrainingPSDs = score(:, 1:numPrincipalComponents);

disp('Post-PCA training set dimensions (samples x components): ');
size(reducedTrainingPSDs)

%% Perform SVM on the learning dataset
disp('Learning SVM model...');

%SVMModel = fitcsvm(comps, classes)
%classOrder = SVMModel.ClassNames

%% Do some predictions on the learning set itself






























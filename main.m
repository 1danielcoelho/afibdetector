%% Cleanup
clc
clear all
close all

%% Import data

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

% Loads all the records into a 'records' structure
for recIndex = 1:length(recordNames)
    recName = char(recordNames(recIndex));
    recPath = strcat('mit-bih\', recName);    
    
    % Actually read everything    
    [r.signalVolts,r.Fs,r.tmSecs]=rdsamp(recPath);
    [r.annSamples,r.anntype,r.subtype,r.chan,r.num,r.comments] = rdann(recPath, 'atr');
    r.annVec = get_af_annotation_vector(length(r.signalVolts), r.annSamples, r.comments);
    
    % Pack everything we read into our records structure
    records.(strcat('rec', recName)) = r;
end

% Go back to where we were before, if it matters
cd(prev_folder);

%% Test plot
testSignal = records.rec04015.signalVolts(103000:110000, 1);
testTime = records.rec04015.tmSecs(103000:110000);

plot(testTime, testSignal);
grid on

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
    disp(recName);
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

%% Test plots
% testSignal = records.rec04015.signalVolts(103000:110000, 1);
% testTime = records.rec04015.tmSecs(103000:110000);
% 
% % Plot some samples of a record
% plot(testTime, testSignal);
% grid on

%% Perform spectrograms
learnRecord = records.rec04015;
testRecord = records.rec04043;

intStart = 103000;
intEnd = 110000;
windowSeconds = 5;  % Chunk size in seconds
overlap = [];  % Samples to overlap, default is 50%
nfft = [];  % How many FFT points desired per chunk. Leaving empty uses window size
samplingFreq = learnRecord.Fs;
[s, f, t, ps] = spectrogram(learnRecord.signalVolts(intStart:intEnd), ...
                     windowSeconds * learnRecord.Fs, ...
                     overlap, ...
                     nfft, ...
                     samplingFreq, ...
                     'yaxis');                 
                 
psdLearn = 10*log10(abs(ps));

[s, f, t, ps] = spectrogram(testRecord.signalVolts(intStart:intEnd), ...
                     windowSeconds * testRecord.Fs, ...
                     overlap, ...
                     nfft, ...
                     samplingFreq, ...
                     'yaxis');                 
                 
psdTest = 10*log10(abs(ps));
                 
%% Test if PSD matches spectrogram's plot
% psd2 = zeros(size(ps, 1)+1, size(ps, 2)+1);
% psd2(1:end-1, 1:end-1) = psdLearn;
% 
% tStep = t(1);
% t2 = 0:tStep:t(end);
% 
% fStep = f(2) * 0.5;
% f2 = 0-fStep:2*fStep:2*fStep*(length(f));
% 
% figure;
% colormap('hot')
% %https://www.mathworks.com/matlabcentral/answers/122472-how-to-get-the-power-spectral-density-from-a-spectrogram-in-a-given-frequency-range
% surf(t2, f2, psd2, 'EdgeColor', 'none')
% colorbar
% view(2)

%% Separate the signal into windows
windowSeconds = 4;
samplingFreq = learnRecord.Fs;
windowSize = windowSeconds * samplingFreq;

% Cada row uma janela
numWindows = floor(length(learnRecord.signalVolts(:, 1)) / windowSize);

classes = zeros(numWindows, 1);
windows = zeros(numWindows, windowSize);

for i=1:numWindows
    rangeStart = 1 + (i-1)*windowSize;
    rangeEnd = i*windowSize;    
    
    windows(i, :) = learnRecord.signalVolts(rangeStart:rangeEnd, 1);    
    sampleAnns = learnRecord.annVec(rangeStart:rangeEnd);
    
    % If most of the samples are marked as AFIB, mark the window as AFIB
    classes(i) = sum(sampleAnns) > windowSize/2;
end

%% Perform PSDs of windows

%https://www.mathworks.com/help/signal/ref/pwelch.html
[psds, f] = pwelch(windows',[],[],[],250);
psds = psds';

% figure;
% plot(f,10*log10(psds(:, 56)));
% xlim([0, 125]);
% title('junto');
% xlabel('Frequency (Hz)')
% ylabel('Magnitude (dB)')

%% Perform LDA of psds
%W = LDA(10*log10(psds), classes);

%% Perform PCA of psds
[coeff,score,latent,tsquared,explained] = pca(10*log10(psds));

% Get first 10 components
comps = score(:, 1:10);






























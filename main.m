%%TODO
% Delay separating training sets from testing sets as much as possible
% Maybe classify between AFIB and normal within the same series (i.e. use
% first 3 hours to train, other 7 to test --> Need to screen for records
% with many episodes)

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

% Each file has two leads
numRecords = length(recordNames) * 2;
           
% Loads all the records into a 'records' structure
for recIndex = 1:length(recordNames)
    recName = char(recordNames(recIndex));
    disp(recName);
    recPath = strcat('mit-bih\', recName);    
    
    % Actually read everything    
    [ecg_signal,r.Fs,r.tmSecs]=rdsamp(recPath);
    [r.annSamples,r.anntype,r.subtype,r.chan,r.num,r.comments] = rdann(recPath, 'atr');
    r.annVec = get_annotation_vector(length(ecg_signal), r.annSamples, r.comments);
        
    % Pack everything we read into our records structure
    r.signalmV = ecg_signal(:, 1);
    records.(strcat('rec', recName, '_1')) = r;    
    
    r.signalmV = ecg_signal(:, 2);
    records.(strcat('rec', recName, '_2')) = r;    
end

% Go back to where we were before, if it matters
cd(prev_folder);

clear ecg_signal

%% Break record signals into windows
disp('Separating ECG1 into windows and extracting their classes...');
windowSizeSeconds = 4;
discardMixedWindows = 1; %Whether we keep windows that are only partially afib/normal

recordNames = fieldnames(records);
for recordName=recordNames'    
    record = records.(recordName{1});
    disp(recordName{1});
    
    samplingFreq = record.Fs;
    PSDSize = windowSizeSeconds * samplingFreq;
        
    maxNumWindows = floor(length(record.signalmV(:, 1)) / PSDSize);
    numPSDs = 1;

    % Pre-allocate maximum size for speed
    labels = zeros(maxNumWindows, 1);
    windows = zeros(maxNumWindows, PSDSize);

    % Break record into windows
    for i=1:maxNumWindows        
        rangeStart = 1 + (i-1)*PSDSize;
        rangeEnd = i*PSDSize;    

        sampleAnns = record.annVec(rangeStart:rangeEnd);
        
        % If there are other arrythmias in this window, discard it and go
        % to the next window
        if sum(sampleAnns==2) > 0
            continue
        end
        
        % If we're asked to discard mixed windows and the number of normal
        % samples is anything other than 0 or the total number of samples,
        % skip this window
        if discardMixedWindows && ((sum(sampleAnns==0) ~= 0) && (sum(sampleAnns==0) ~= length(sampleAnns)))
            continue
        end
        
        % If most of the samples are marked as AFIB, mark the window as AFIB        
        labels(numPSDs) = sum(sampleAnns) > PSDSize/2;
        windows(numPSDs, :) = record.signalmV(rangeStart:rangeEnd, 1);        
        
        numPSDs = numPSDs + 1;
    end    
    
    % Rewind last iteration just before we left the for loop
    numPSDs = numPSDs - 1;
    
    % Discard extra lines
    windows = windows(1:numPSDs, :);
    labels = labels(1:numPSDs);
    
    records.(recordName{1}).signalmVWindows = windows;
    records.(recordName{1}).actualClasses = labels;
    
    % Delete variables we won't use anymore
    records.(recordName{1}) = rmfield(records.(recordName{1}), 'signalmV');
    records.(recordName{1}) = rmfield(records.(recordName{1}), 'annVec');
end

clear windows

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
    
    % Check for Infs and NaNs (shouldn't happen anymore)
    nans = sum(isnan(records.(recordName{1}).PSDs));
    infs = sum(isinf(records.(recordName{1}).PSDs));    
    if nans
        disp(strcat(recordName{1}, ' has NaNs in its PSDs'));
    end
    
    if infs
        disp(strcat(recordName{1}, ' has Infs in its PSDs'));
    end
end

%% Discard records that don't have enough windows with normal/afib waveforms
disp('Discarding records with too few windows of normal or afib... ');
minNumWindowsPerRecord = 50;

recordNames = fieldnames(records);
for recordName=recordNames'     
    record = records.(recordName{1});    
    
    normalWindows = sum(record.actualClasses == 0);
    afibWindows = sum(record.actualClasses == 1);
    
    % Print stuff for a table later
    %disp(strcat(recordName{1}, ',', int2str(normalWindows), ',', int2str(afibWindows)));
    
    if normalWindows < minNumWindowsPerRecord || afibWindows < minNumWindowsPerRecord
        disp(recordName{1});
        records = rmfield(records, recordName{1});
    end
end

%% Get a sample of a normal and afib window from each record, as well as their PSDs
recordNames = fieldnames(records);

normalWindowPSD = zeros(length(recordNames), size(psds, 2));  
AFIBWindowPSD = zeros(length(recordNames), size(psds, 2));  
recordIndex = 1;

for recordName=recordNames'     
    record = records.(recordName{1});    
    
    timeAxis = 0:(1.0/record.Fs):windowSizeSeconds - 1/record.Fs;
    
    windows = record.signalmVWindows;
    labels = record.actualClasses;
    psds = record.PSDs;
        
    normalWindowIndices = find(~labels);
    AFIBWindowIndices = find(labels);
    
    % Get a normal and an AFIB window to plot as examples
    firstNormalWindow = windows(normalWindowIndices(1), :);
    firstAFIBWindow = windows(AFIBWindowIndices(1), :);
    
    % Get mean of all normal window PSDs
    for i=1:length(normalWindowIndices)
        normalWindowPSD(recordIndex, :) = normalWindowPSD(recordIndex, :) + psds(normalWindowIndices(i), :);
    end
    normalWindowPSD(recordIndex, :) = normalWindowPSD(recordIndex, :)./length(normalWindowIndices);    
    
    % Get mean of all AFIB window PSDs
    for i=1:length(AFIBWindowIndices)
        AFIBWindowPSD(recordIndex, :) = AFIBWindowPSD(recordIndex, :) + psds(AFIBWindowIndices(i), :);
    end
    AFIBWindowPSD(recordIndex, :) = AFIBWindowPSD(recordIndex, :)./length(AFIBWindowIndices);    
    
    % Plot everything
    figure;
    suptitle(strcat(recordName{1}));
    
    subplot(2, 2, 1);
    plot(timeAxis, firstNormalWindow);
    title('First normal window');
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');
    
    subplot(2, 2, 2);
    plot(timeAxis, firstAFIBWindow);
    title('First AFIB window');
    xlabel('Time (s)');
    ylabel('Amplitude (mV)');    
    
    subplot(2, 2, 3:4);    
    plot(normalWindowPSD(recordIndex, :));
    hold on
    plot(AFIBWindowPSD(recordIndex, :));    
    title('Average of window PSDs');
    xlabel('Frequency (Hz)');
    ylabel('Power Spectral Density (dB)');    
    legend('Normal', 'AFIB');
end

% Average all normal PSDs and all AFIB PSDs from all records
figure;
suptitle('Mean PSD of all windows from all records');  
plot(mean(normalWindowPSD));
hold on
plot(mean(AFIBWindowPSD));    
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB)');    
legend('Normal', 'AFIB');

clear timeAxis
clear normalWindowIndices AFIBWindowIndices
clear psds windows labels
clear normalWindowPSD AFIBWindowPSD
clear firstNormalWindow firstAFIBWindow

%% Delete signal and time, since we won't use them anymore
recordNames = fieldnames(records);
for recordName=recordNames'    
    records.(recordName{1}) = rmfield(records.(recordName{1}), 'signalmVWindows');
    records.(recordName{1}) = rmfield(records.(recordName{1}), 'tmSecs');
end

%% Extract training windows from the data
disp('Separating training windows...');

% We will use this percentage of the psds as a training set
trainingPercentage = 0.3;

recordNames = fieldnames(records);
for recordName=recordNames'     
    record = records.(recordName{1});
    
    labels = logical(record.actualClasses);
    numPSDs = size(labels, 1);
    PSDSize = size(record.PSDs, 2);    
    
    % Determine how many normal/AFIB training/test windows we'll need
    numNormalPSDs = sum(labels==0);
    numTrainingNormalPSDs = min(max(round(trainingPercentage * numNormalPSDs), 1), numNormalPSDs-1);
    numTestNormalPSDs = numNormalPSDs - numTrainingNormalPSDs;        
    numAFIBPSDs = sum(labels==1);            
    numTrainingAFIBPSDs = min(max(round(trainingPercentage * numAFIBPSDs), 1), numAFIBPSDs-1);
    numTestAFIBPSDs = numAFIBPSDs - numTrainingAFIBPSDs;
    
    % Separate all normal windows from AFIB windows
    normalPSDs = record.PSDs(~labels, :);
    AFIBPSDs = record.PSDs(labels, :);    
    
    % Shuffle normal and AFIB PSDs so we pick random samples to be training
    % sets and testing sets    
    normalPSDs = normalPSDs(randperm(numNormalPSDs), :);
    AFIBPSDs = AFIBPSDs(randperm(numAFIBPSDs), :);
    
    % Get the four groups of PSDs 
    trainingNormalPSDs = normalPSDs(1:numTrainingNormalPSDs, :);
    testNormalPSDs = normalPSDs(numTrainingNormalPSDs+1:end, :);    
    trainingAFIBPSDs = AFIBPSDs(1:numTrainingAFIBPSDs, :);    
    testAFIBPSDs = AFIBPSDs(numTrainingAFIBPSDs+1:end, :);
    
    % Combine them into training and testing sets
    trainingPSDs = [trainingNormalPSDs; trainingAFIBPSDs];
    testPSDs = [testNormalPSDs; testAFIBPSDs];
    
    % Create training and test class labels
    trainingClasses = true(numTrainingNormalPSDs + numTrainingAFIBPSDs, 1);
    trainingClasses(1:numTrainingNormalPSDs) = 0;    
    testClasses = true(numTestNormalPSDs + numTestAFIBPSDs, 1);
    testClasses(1:numTestNormalPSDs) = 0;
    
    % Pack them back into the record
    records.(recordName{1}).trainingPSDs = trainingPSDs;
    records.(recordName{1}).testPSDs = testPSDs;
    records.(recordName{1}).trainingClasses = trainingClasses;
    records.(recordName{1}).testClasses = testClasses;       
    
    % Delete variables we won't use anymore
    records.(recordName{1}) = rmfield(records.(recordName{1}), 'PSDs');
    records.(recordName{1}) = rmfield(records.(recordName{1}), 'actualClasses');
end

clear AFIBPSDs normalPSDs
clear testAFIBPSDs testNormalPSDs
clear trainingAFIBPSDs trainingNormalPSDs
clear trainingClasses testClasses
clear trainingPSDs testPSDs 

%% Perform PCA of training PSDs
% 
% numPrincipalComponents = 0;
% 
% % Only perform PCA if we pick numPrincipalComponents different than zero
% if numPrincipalComponents ~= 0
%     disp('Performing PCA from window PSDs...');
%     
%     [coeff,score,latent,tsquared,explained] = pca(trainingPSDs);
% 
%     % We'll use this later to convert new observations into PCA components
%     pcaMean = mean(trainingPSDs);
%     pcaCoeffs = inv(coeff');
% 
%     % Get first 10 components
%     trainingPSDsPCA = score(:, 1:numPrincipalComponents);
% 
%     % Map test records to PCA components
%     recordNames = fieldnames(records);
%     for recordName=recordNames'     
%         record = records.(recordName{1});
%         testPSDs = record.PSDs;        
%         numSamples = size(testPSDs, 1);
% 
%         if ~record.isLearningSet       
%             % Actually do the mapping of samples. Check test_pca for proof
%             meanPCAmat = repmat(pcaMean, numSamples, 1);
%             pcaTestPSDs = (testPSDs - meanPCAmat) * pcaCoeffs;
%             pcaTestPSDs = pcaTestPSDs(:, 1:numPrincipalComponents);
%             
%             records.(recordName{1}).PSDs = pcaTestPSDs;
%         end
%     end
% 
%     disp('Post-PCA training set dimensions (samples x components): ');
%     size(trainingPSDsPCA)
%     
%     % Overwrite our training PSDs so that the rest of the pipeline never
%     % needs to care if we did PCA or not
%     trainingPSDs = trainingPSDsPCA;
% end

%% Try a linear classifier instead
disp('Classifying with a linear classifier...');

totalConfMat = zeros(2, 2);

recordNames = fieldnames(records);
for recordName=recordNames'     
    record = records.(recordName{1});    
      
%     % Linear classifier
%     predictedClasses = classify(record.testPSDs, ...
%                                 record.trainingPSDs, ...
%                                 record.trainingClasses, ...
%                                 'linear');

    % SVM Classifier
    SVMModel = fitcsvm(record.trainingPSDs, ...
                       record.trainingClasses, ...
                       'KernelFunction', ...
                       'polynomial', ...
                       'PolynomialOrder', ...
                       2);
    [predictedClasses, scoreForEachClass] = predict(SVMModel, record.testPSDs);
                   
    actualClasses = record.testClasses;

    confMat = zeros(2, 2);        
    confMat(1, 1) = sum(~predictedClasses & ~actualClasses);
    confMat(2, 2) = sum(predictedClasses & actualClasses);
    confMat(1, 2) = sum(predictedClasses & ~actualClasses);
    confMat(2, 1) = sum(~predictedClasses & actualClasses);        

%     disp(recordName{1});
%     accuracy = (confMat(1, 1) + confMat(2, 2)) / sum(sum(confMat))
    
    totalConfMat = totalConfMat + confMat;

    % Store prediction results into records
    records.(recordName{1}).predictedClass = predictedClasses;
end

totalConfMat
totalAccuracy = (totalConfMat(1, 1) + totalConfMat(2, 2)) / sum(sum(totalConfMat))























function [ psds ] = get_psd_from_record( record, windowSeconds )
%get_psd_from_record Returns PSD according to Welch estimator of the
%signalVolts contained in record

samplingFreq = record.Fs;
windowSize = windowSeconds * samplingFreq;
numWindows = floor(length(record.signalVolts(:, 1)) / windowSize);

classes = zeros(numWindows, 1);
windows = zeros(numWindows, windowSize);

% Break record into windows
for i=1:numWindows
    rangeStart = 1 + (i-1)*windowSize;
    rangeEnd = i*windowSize;    
    
    windows(i, :) = record.signalVolts(rangeStart:rangeEnd, 1);    
    sampleAnns = record.annVec(rangeStart:rangeEnd);
    
    % If most of the samples are marked as AFIB, mark the window as AFIB
    classes(i) = sum(sampleAnns) > windowSize/2;
end

% Obtain PSD with Welch estimator https://www.mathworks.com/help/signal/ref/pwelch.html
[psds, f] = pwelch(windows',[],[],[],250);
psds = psds';

end


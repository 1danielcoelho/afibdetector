%Source: https://www.mathworks.com/matlabcentral/answers/53259-how-to-project-a-new-point-to-pca-new-basis

load hald;
data = ingredients;

% PCA, W is coefficient and Y is the score
[coef,score] = pca(data);

%% Map a group of observations to PCA component space manually
% According to the reconstruction rule, we should have x1=y1*W'
% therefore, y1 = x1/W'

% Get a group of entries
g = data(1:4, :);

% Batch convert
numSamples = size(g, 1);
pcaMean = mean(data);
pcaCoeffs = inv(coef');

meanPCAmat = repmat(pcaMean, numSamples, 1);
y = (g - meanPCAmat) * pcaCoeffs

% Compare with the converted values that PCA itself returns
score(1:4, :)









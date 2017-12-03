%https://www.mathworks.com/matlabcentral/answers/53259-how-to-project-a-new-point-to-pca-new-basis

load hald;
data = ingredients;

% Weight and Mean:

% PCA, W is coefficient and Y is the score
[coef,score] = pca(data);

% First observation of the centered data and its score
x1 = data(1,:) - mean(data);
y1 = score(1,:)

% According to the reconstruction rule, we should have x1=y1*W'
% therefore, y1 = x1/W'
y = x1/coef'

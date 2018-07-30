
% Add the ImageGraph to path (in a folder named 'ImageGraphs'). Find it here:
% http://www.mathworks.com/matlabcentral/fileexchange/53614-image-graphs
addpath(fullfile(fileparts(mfilename('fullpath')), 'ImageGraphs'));

%% Read the image and features from the sample file
image = im2double(imread('docia.png'));
features = image(:, size(image, 2) / 2 + 1 : end, :);
image = image(:, 1 : size(image, 2) / 2, :);

% The eigendecomposition uses a lot of memory and may render the computer
% unresponsive, so better to test it first with a small image.
image = imresize(image, 0.5);
features = imresize(features, 0.5);

%% Semantic soft segmentation
% This function outputs many intermediate variables, if needed.
% The results may vary a bit from run to run, as there are 2 stages that use 
% k-means for intialization & grouping.
sss = SemanticSoftSegmentation(image, features);

% To use the features generated using our network implementation,
% just feed them as the 'features' variable to the function. It will do
% the prepocessing described in the paper and give the processed
% features as an output.
% If you are dealing with many images, storing the features after
% preprocessing is recommended as raw hyperdimensional features
% take a lot of space. Check the 'preprocessFeatures.m' file.

% Visualize
figure; imshow([image features visualizeSoftSegments(sss)]);
title('Semantic soft segments');

% There's also an implementation of Spectral Matting included
sm = SpectralMatting(image);
% You can group the soft segments from Spectral Matting using
% semantic features, the way we presented our comparisons in the paper.
sm_gr = groupSegments(sm, features);
figure; imshow([image visualizeSoftSegments(sm) visualizeSoftSegments(sm_gr)]);
title('Matting components');

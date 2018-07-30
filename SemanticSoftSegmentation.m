
% Semantic Soft Segmentation
% This function implements the soft segmentation approach described in
% Yagiz Aksoy, Tae-Hyun Oh, Sylvain Paris, Marc Pollefeys, Wojciech Matusik
% "Semantic Soft Segmentation", ACM TOG (Proc. SIGGRAPH) 2018

function [softSegments, initSoftSegments, Laplacian, affinities, features, superpixels, eigenvectors, eigenvalues] = SemanticSoftSegmentation(image, features)
    
    disp('Semantic Soft Segmentation')
    % Prepare the inputs and superpixels
    image = im2double(image);
    if size(features, 3) > 3 % If the features are raw, hyperdimensional, preprocess them
        features = preprocessFeatures(features, image);
    else
        features = im2double(features);
    end
    superpixels = Superpixels(image);
    [h, w, ~] = size(image);

    disp('     Computing affinities')
    % Compute the affinities and the Laplacian
    affinities{1} = mattingAffinity(image);
    affinities{2} = superpixels.neighborAffinities(features); % semantic affinity
    affinities{3} = superpixels.nearbyAffinities(image); % non-local color affinity
    Laplacian = affinityMatrixToLaplacian(affinities{1} + 0.01 * affinities{2} + 0.01 * affinities{3}); % Equation 6

    disp('     Computing eigenvectors')
    % Compute the eigendecomposition
    eigCnt = 100; % We use 100 eigenvectors in the optimization
    [eigenvectors, eigenvalues] = eigs(Laplacian, eigCnt, 'SM');
    
    disp('     Initial optimization')
    % Compute initial soft segments
    initialSegmCnt = 40;
    sparsityParam = 0.8;
    iterCnt = 40;
    % feeding features to the function below triggers semantic intialization
    initSoftSegments = softSegmentsFromEigs(eigenvectors, eigenvalues, Laplacian, ...
                                            h, w, features, initialSegmCnt, iterCnt, sparsityParam, [], []);

    % Group segments w.r.t. their mean semantic feature vectors
    groupedSegments = groupSegments(initSoftSegments, features);

    disp('     Final optimization')
    % Do the final sparsification
    softSegments = sparsifySegments(groupedSegments, Laplacian, imageGradient(image, false, 6));
    
    disp('     Done.')
end
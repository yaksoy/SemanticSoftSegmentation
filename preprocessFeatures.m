
% Pre-processing hyper-dimensional semantic feature vectors
% as described in Section 3.5.

function simp = preprocessFeatures(features, image)
    % Filter out super high numbers due to some instability in the network
    features(features < -5) = -5;
    features(features > 5) = 5;
    
    % Filter each channel of features with image as the guide
    if exist('image', 'var') && ~isempty(image)
        fd = size(features, 3);
        maxfd = fd - rem(fd, 3);
        for i = 1 : 3 : maxfd
            features(:, :, i : i+2) = imguidedfilter(features(:, :, i : i+2), image, 'NeighborhoodSize', 10);
        end
        for i = maxfd + 1 : fd
            features(:, :, i) = imguidedfilter(features(:, :, i), image, 'NeighborhoodSize', 10);
        end
    end

    % Run PCA and normalize to [0, 1]
    simp = featuresPCA(features, 3);
    for i = 1 : 3
        simp(:,:,i) = simp(:,:,i) - min(min(simp(:,:,i)));
        simp(:,:,i) = simp(:,:,i) / max(max(simp(:,:,i)));
    end
end

function pcafeat = featuresPCA(features, dim)
    features = double(features);
    [h, w, d] = size(features);
    features = reshape(features, [h*w, d]);
    featmean = mean(features, 1);
    features = features - ones(h*w, 1) * featmean;
    covar = features' * features;
    [eigvecs, ~] = eigs(covar, dim, 'LA');
    pcafeat = features * eigvecs;
    pcafeat = reshape(pcafeat, [h, w, dim]);
end

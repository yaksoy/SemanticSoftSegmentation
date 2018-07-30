
% Spectral Matting
% This function implements the soft segmentation approach described in
% Anat Levin, Dani Lischinski, Yair Weiss, "Spectral Matting", IEEE TPAMI, 2008.
% The parameters here are set to their default values as reported by Levin et al.

function [softSegments, Laplacian, eigenvectors, eigenvalues] = SpectralMatting(image)
    
    disp('Spectral Matting')
    image = im2double(image);
    [h, w, ~] = size(image);

    disp('     Computing affinities')
    % Compute the matting Laplacian
    Laplacian = affinityMatrixToLaplacian(mattingAffinity(image));

    disp('     Computing eigenvectors')
    % Compute the eigendecomposition
    eigCnt = 50;
    [eigenvectors, eigenvalues] = eigs(Laplacian, eigCnt, 'SM');
    
    disp('     Optimization')
    % Compute the soft segments = matting components
    initialSegmCnt = 20;
    sparsityParam = 0.8;
    iterCnt = 20;
    softSegments = softSegmentsFromEigs(eigenvectors, eigenvalues, Laplacian, ...
                                            h, w, [], initialSegmCnt, iterCnt, sparsityParam, [], []);
                                            
    disp('     Done.')
end

function vis = visualizeEigenvectorRedGreen(eigenvector)
    % negative values are shown in red, and positive in green
    vis = zeros(size(eigenvector, 1), size(eigenvector, 2), 3);
    vis(:,:,1) = -100 * eigenvector;
    vis(:,:,2) = 100 * eigenvector;
    vis(vis < 0) = 0;
    vis(vis > 1) = 1;
end

% This function implements the relaxed sparsification descibed in Section 3.4

function softSegments = sparsifySegments(softSegments, Laplacian, imageGrad)
 
    sigmaS = 1; % sparsity
    sigmaF = 1; % fidelity
    delta = 100; % constraint
    [h, w, compCnt] = size(softSegments);
    N = h * w * compCnt;

    if ~exist('imageGrad', 'var') || isempty(imageGrad)
        % If image gradient is not provided, set the param to the default 0.9
        spPow = 0.90;
    else
        % Compute the per-pixel sparsity parameter from the gradient
        imageGrad(imageGrad > 0.1) = 0.1;
        imageGrad = imageGrad + 0.9;
        spPow = repmat(imageGrad(:), [compCnt, 1]);
    end

    % Iter count for pcg and main optimization
    itersBetweenUpdate = 100;
    highLevelIters = 20;

    % Get rid of very low/high alpha values and normalize
    softSegments(softSegments < 0.1) = 0;
    softSegments(softSegments > 0.9) = 1;
    softSegments = softSegments ./ repmat(sum(softSegments, 3), [1 1 size(softSegments, 3)]);
    
    % Construct the linear system
    lap = Laplacian;
    for i = 2 : compCnt
        Laplacian = blkdiag(Laplacian, lap);
    end

    % The alpha constraint
    C = repmat(speye(h*w), [1 compCnt]);
    C = C' * C;
    Laplacian = Laplacian + delta * C;

    % The sparsification optimization
    softSegments = softSegments(:);
    compInit = softSegments; % needed for fidelity energy
    for iter = 1 : highLevelIters
        if rem(iter, 5) == 0
            disp(['               Iteration ' int2str(iter) ' of ' int2str(highLevelIters)]);
        end
        [u, v] = getUandV(softSegments, spPow); % The sparsity energy
        A = Laplacian + sigmaS * (spdiags(u, 0, N, N) + spdiags(v, 0, N, N)) + sigmaF * speye(N);
        b = sigmaS * v + sigmaF * compInit + delta;
        [softSegments, ~] = pcg(A, b, [], itersBetweenUpdate, [], [], softSegments);
    end

    % One final iter for good times (everything but sparsity)
    A = Laplacian + sigmaF * speye(N);
    b = sigmaF * softSegments + delta;
    softSegments = pcg(A, b, [], 10 * itersBetweenUpdate, [], [], softSegments);

    % Ta-dah
    softSegments = reshape(softSegments, [h w compCnt]);
end

function [u, v] = getUandV(comp, spPow)
    % Sparsity terms in the energy
    eps = 1e-2;
    u = max(abs(comp(:)), eps) .^ (spPow - 2);
    v = max(abs(1 - comp(:)), eps) .^ (spPow - 2);
end
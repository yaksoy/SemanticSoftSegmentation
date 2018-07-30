
function Lap = affinityMatrixToLaplacian(aff, normalize)
    if ~exist('normalize', 'var') || isempty(normalize)
        normalize = false ;
    end
    N = size(aff, 1);
    if normalize
        aa = sum(aff, 2);
        D = spdiags(aa, 0 , N, N);
        aa = sqrt(1./aa);
        D12 = spdiags(aa, 0 , N, N);
        Lap = D12 * (D - aff) * D12;
    else
        Lap = spdiags(sum(aff, 2), 0 , N, N) - aff;
    end
end
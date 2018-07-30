
% This function implements the constrained sparsification descibed in Section 3.4.
% This methodology is intriduced by Levin et al. in "Spectral Matting", and
% this function is an adapted version of their original source code. Find it here:
% http://www.vision.huji.ac.il/SpectralMatting/

function [softSegments, initialSegments] = softSegmentsFromEigs(eigVecs, eigVals, Laplacian, h, w, features, compCnt, maxIter, sparsityParam, imageGrad, initialSegments)
    if (~exist('maxIter','var') || isempty(maxIter))
      maxIter = 20;
    end
    if (~exist('sparsityParam','var') || isempty(sparsityParam))
        sparsityParam = 0.9;
    end
    if (~exist('compCnt','var') || isempty(compCnt))
        compCnt = 10;
    end
    eigValCnt = size(eigVecs, 2);
    if (~exist('eigVals','var') || isempty(eigVals))
        eigVals = 1e-5*eye(eigValCnt);
    end

    % Initialize with k-means if an initialization (initialSegments) isn't provided
    if ~exist('initialSegments', 'var') || isempty(initialSegments)
      % Use the features for the k-means, if provided
      if ~(~exist('features', 'var') || isempty(features))
        disp('          Computing k-means initialization using semantic features...');
        %features = featuresPCA(features, 10);
        features = reshape(features, [size(features, 1) * size(features, 2), size(features, 3)]);
        [initialSegments, ~, ~, ~] = fastKmeans(features, compCnt);
      else
      % This is the default initialization from Spectral Matting
        disp('          Computing k-means initialization...');
        eigVals = abs(eigVals); % Sometimes there are -epsilon numbers that mess up the thing
        initEigsCnt = 20;
        initEigsWeights = diag(1 ./ diag(eigVals(2 : initEigsCnt + 1, 2 : initEigsCnt + 1) .^ 0.5 ));
        initEigVecs = eigVecs(:, 2 : initEigsCnt + 1) * initEigsWeights;
        [initialSegments, ~, ~, ~] = fastKmeans(initEigVecs, compCnt);
      end
    end
    softSegments = zeros(length(initialSegments), compCnt);
    for i = 1 : compCnt
        softSegments(:, i) = double(initialSegments == i);
    end

    if ~exist('imageGrad', 'var') || isempty(imageGrad)
        spMat = sparsityParam;
    else
        imageGrad(imageGrad > 0.2) = 0.2;
        imageGrad = imageGrad + 0.8;
        spMat = repmat(imageGrad(:), [1, compCnt]);
    end

    % Second derivative of first and second terms in the sparsity penalty
    disp('          Starting optimization...');
    thr_e = 1e-10;
    w1 = 0.3;
    w0 = 0.3;
    e1 = w1 .^ sparsityParam * max(abs(softSegments-1), thr_e) .^ (spMat - 2);
    e0 = w0 .^ sparsityParam * max(abs(softSegments), thr_e) .^ (spMat - 2);

    scld = 1;
    eig_vectors = eigVecs(:, 1 : eigValCnt);
    eig_values = eigVals(1 : eigValCnt, 1 : eigValCnt);

    % First iter no for removing zero components
    removeIter = ceil(maxIter / 4);
    removeIterCycle = ceil(maxIter / 4);

    % Compute matting component with sparsity prior
    for iter = 1 : maxIter
        if rem(iter, 10) == 0
            disp(['               Iteration ' int2str(iter) ' of ' int2str(maxIter)]);
        end
        % Construct the matrices in Eq(9) in Spectral Matting
        tA = zeros((compCnt - 1) * eigValCnt);
        tb = zeros((compCnt - 1) * eigValCnt, 1);
        for k = 1 : compCnt - 1
            weighted_eigs = repmat(e1(:, k) + e0(:, k), 1, eigValCnt) .* eig_vectors;
            tA((k-1) * eigValCnt + 1 : k * eigValCnt, (k-1) * eigValCnt + 1 : k * eigValCnt) = eig_vectors' * weighted_eigs + scld * eig_values;
            tb((k-1) * eigValCnt + 1 : k * eigValCnt) = eig_vectors' * e1(:,k);
        end 
        k = compCnt;
        weighted_eigs = repmat(e1(:, k) + e0(:, k), 1, eigValCnt) .* eig_vectors;
        ttA = eig_vectors' * weighted_eigs + scld * eig_values;
        ttb = eig_vectors' * e0(:, k) + scld * sum(eig_vectors' * Laplacian, 2);

        tA = tA + repmat(ttA, [compCnt - 1, compCnt - 1]);
        tb = tb + repmat(ttb, [compCnt - 1, 1]);

        % Solve for weights
        y = reshape(tA \ tb, eigValCnt, compCnt - 1);

        % Compute the matting comps from weights
        softSegments = eigVecs(:, 1 : eigValCnt) * y;
        softSegments(:, compCnt) = 1 - sum(softSegments(:, 1 : compCnt - 1), 2); % Sets the last one as 1-sum(others), guaranteeing \sum(all) = 1

        % Remove matting components which are close to zero every once in a while
       if iter > removeIter
           nzii = find(max(abs(softSegments)) > 0.1);
           compCnt = length(nzii);
          softSegments = softSegments(:, nzii);
           removeIter = removeIter + removeIterCycle;
       end

        % Recompute the derivatives of sparsity penalties
        if length(spMat(:)) == 1
          e1 = w1 .^ sparsityParam * max(abs(softSegments-1), thr_e) .^ (spMat - 2);
          e0 = w0 .^ sparsityParam * max(abs(softSegments), thr_e) .^ (spMat - 2);
        else
          e1 = w1 .^ sparsityParam * max(abs(softSegments-1), thr_e) .^ (spMat(:, 1 : size(softSegments, 2)) - 2);
          e0 = w0 .^ sparsityParam * max(abs(softSegments), thr_e) .^ (spMat(:, 1 : size(softSegments, 2)) - 2);
        end
    end
    softSegments = reshape(softSegments, [h, w, size(softSegments, 2)]);
end

function [idx, C, sumd, D]=fastKmeans(X,K)
  % X: points in the N-by-P data matrix
  % idx - an N-by-1 vector containing the cluster indices of each point
  %
  % c - the K cluster centroid locations in a K-by-P matrix.
  % sumd - the within-cluster sums of point-to-centroid distances in a 1-by-K vector.
  % distances from each point to every centroid in a N-by-K.

  startK = 5;
  startK = min(K,startK);
  maxIters = 100;  % Defualt of matlab is 100    

    X = sign(real(X)) .* abs(X);
  
  [idx, C, sumd, D]=kmeans(X,startK,'EmptyAction','singleton','Start','cluster', ...
      'Maxiter', maxIters, 'Replicates', 7);

  valid_vec = zeros(1,startK);
  scr_vec = zeros(1,startK)-1;

  for compCnt = startK+1:K
      % create a new cluster by splitting each cluster to two...
      max_scr=-1;
      clear min_C;
      for cl = 1:compCnt-1
        cl_mask = idx == cl;
        cl_idxs = find(cl_mask);
        clX = X(cl_idxs,:);
        if (size(clX,1)> 2*size(clX,2))
          if (valid_vec(cl) == 0)
            [tmp_idx, tmp_C, ~, tmp_D]=kmeans(clX,2,'EmptyAction','singleton','Start','cluster', 'Maxiter', maxIters);
            % chk how much the partition helps ...
            scr=sum(min(D(cl_idxs,:),[],2))-sum(min(tmp_D,[],2));
            scr_vec(cl) = scr;
          else % we already saved it...
            scr = scr_vec(cl);
          end         
        else
          scr=-2;
          scr_vec(cl) = scr;
        end  
        if (scr > max_scr)
          if (valid_vec(cl)==1) % not for the scr. Just for the idxs.
            [tmp_idx, tmp_C, ~, ~]=kmeans(clX,2,'EmptyAction','singleton','Start','cluster', 'Maxiter', maxIters);
          end
          
          max_scr = scr;
          bestC = [C;tmp_C(2,:)];  
          bestC(cl,:) = tmp_C(1,:);
          best_cl = cl;         

          best_idx = idx;
          best_idx(cl_idxs) = (tmp_idx == 1)*best_cl + (tmp_idx == 2)*compCnt;
        end
        valid_vec(cl) = 1;
      end        
      C = bestC;
      idx = best_idx;
      
      valid_vec = [valid_vec, 0];   % the two new clusers are new, so their
      valid_vec(best_cl) = 0;      % score have not been computed yet.
      scr_vec = [scr_vec, -1];
      scr_vec(best_cl) = -1;

      if (compCnt < 13)
        [idx, C, sumd, D]=kmeans(X, compCnt, 'EmptyAction', 'singleton', 'Start', C, 'Maxiter', maxIters);       
        valid_vec = zeros(1,compCnt);
      end
  end

end
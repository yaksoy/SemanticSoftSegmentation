
% This class is for computing superpixel-based affinities described in the paper.
% This class requires the image graphs methods by Steve Eddins. Find it here:
% http://www.mathworks.com/matlabcentral/fileexchange/53614-image-graphs

classdef Superpixels
    properties 
        labels
        spcount
        neigh
        centroids
    end
    methods
        function obj = Superpixels(im, spcnt)
            if ~exist('spcnt', 'var') || isempty(spcnt)
                spcnt = 2500;
            end
            [L, N] = superpixels(im, spcnt, 'Compactness', 1e-20);
            obj.labels = L;
            obj.spcount = N;
            % Find neighboring superpixels
            g = adjacentRegionsGraph(L);
            obj.neigh = g.Edges.Labels;
            % Find centroids
            s = regionprops(L, 'centroid');
            cent = cat(1, s.Centroid);
            obj.centroids = round(cent(:, 2:-1:1));
            [h, w, ~] = size(im);
            obj.centroids(:, 3) = sub2ind([h, w], obj.centroids(:, 1), obj.centroids(:, 2));
        end

        function regmeans = computeRegionMeans(obj, image)
            [h, w, c] = size(image);
            image = reshape(image, [h*w, c]);
            regmeans = zeros(obj.spcount, c);
            idx = label2idx(obj.labels);
            for i = 1 : length(idx)
                regmeans(i, :) = mean(image(idx{i}, :), 1);
            end
        end

        % This is for the semantic affinity, generates affinities in [-1, 1]
        function W = neighborAffinities(obj, features, erfSteepness, erfCenter)
            if ~exist('erfSteepness', 'var') || isempty(erfSteepness)
                erfSteepness = 20;
            end
            if ~exist('erfCenter', 'var') || isempty(erfCenter)
                erfCenter = 0.85;
            end
            [h, w, ~] = size(features);
            N = h * w;
            spMeans = obj.computeRegionMeans(features);
            affs = zeros(size(obj.neigh, 1), 1);
            inds1 = affs;
            inds2 = affs;
            for i = 1 : size(obj.neigh, 1)
                ind1 = obj.neigh(i, 1);
                ind2 = obj.neigh(i, 2);
                affs(i) = sigmoidAff(spMeans(ind1, :), spMeans(ind2, :), erfSteepness, erfCenter);
                inds1(i) = obj.centroids(ind1, 3);
                inds2(i) = obj.centroids(ind2, 3);
            end
            W = sparse(inds1, inds2, affs, N, N);
            W = W' + W;
        end
        
        % This is for the nonlocal color affinity, generates affinities in [0, 1]
        function W = nearbyAffinities(obj, image, erfSteepness, erfCenter, proxThresh)
            if ~exist('erfSteepness', 'var') || isempty(erfSteepness)
                erfSteepness = 50;
            end
            if ~exist('erfCenter', 'var') || isempty(erfCenter)
                erfCenter = 0.95;
            end
            if ~exist('proxThresh', 'var') || isempty(proxThresh)
                proxThresh = 0.2;
            end
            [h, w, ~] = size(image);
            N = h * w;
            spMeans = obj.computeRegionMeans(image);
            combinationCnt = obj.spcount;
            combinationCnt = combinationCnt * (combinationCnt - 1) / 2;
            affs = zeros(combinationCnt, 1);
            inds1 = affs;
            inds2 = affs;
            cnt = 1;
            cents = obj.centroids(:, 1:2);
            cents(:,1) = cents(:,1) / h;
            cents(:,2) = cents(:,2) / w;
            for i = 1 : obj.spcount
                for j = i + 1 : obj.spcount
                    centdist = cents(i, 1:2) - cents(j, 1:2);
                    centdist = sqrt(centdist * centdist');
                    if centdist > proxThresh
                        affs(cnt) = 0;
                    else
                        affs(cnt) = sigmoidAffPos(spMeans(i, :), spMeans(j, :), erfSteepness, erfCenter);
                    end
                    inds1(cnt) = obj.centroids(i, 3);
                    inds2(cnt) = obj.centroids(j, 3);
                    cnt = cnt + 1;
                end
            end
            W = sparse(inds1, inds2, affs, N, N);
            W = W' + W;
        end

        function vis = visualizeRegionMeans(obj, im)
            vis = label2rgb(obj.labels, obj.computeRegionMeans(im));
        end

    end
end

function aff = sigmoidAff(feat1, feat2, steepness, center)
    aff = abs(feat1 - feat2);
    aff = 1 - sqrt(aff * aff');
    aff = (erf(steepness * (aff - center)));
end

function aff = sigmoidAffPos(feat1, feat2, steepness, center)
    aff = abs(feat1 - feat2);
    aff = 1 - sqrt(aff * aff');
    aff = (erf(steepness * (aff - center)) + 1) / 2;
end
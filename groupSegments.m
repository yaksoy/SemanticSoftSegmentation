
% A simple grouping of soft segments w.r.t. their semantic features
% as described in Section 3.4.

function groupedSegments = groupSegments(segments, features, segmCnt)
    if ~exist('segmCnt', 'var') || isempty(segmCnt)
        segmCnt = 5;
    end
    [h, w, cnt] = size(segments);
    compFeatures = zeros(cnt, size(features, 3));
    for i = 1 : cnt
        cc = segments(:,:,i) .* features;
        cc = sum(sum(cc, 1), 2) / sum(sum(segments(:,:,i), 1), 2);
        compFeatures(i, :) = cc;
    end
    ids = kmeans(compFeatures, segmCnt);
    groupedSegments = zeros(h, w, segmCnt);
    for i = 1 : segmCnt
        groupedSegments(:,:,i) = sum(segments(:,:,ids==i), 3);
    end
end
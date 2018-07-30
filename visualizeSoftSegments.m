
% Assigns a distinct solid color to eact soft segment and composites
% them using the corresponding alpha calues

function [vis, softSegments] = visualizeSoftSegments(softSegments, doOrdering)
    if ~exist('doOrdering', 'var') || isempty(doOrdering)
        doOrdering = false;
    end
    if doOrdering
        % Order layers w.r.t. sum(alpha(:)) -- makes visualizations more consistent across images
        order = determineAlphaOrder(softSegments);
        softSegments = softSegments(:, :, order);
    end
    % A hard-coded set of 'distinct' colors
    colors = createPalette();
    % One solid color per segment, mixed w.r.t. alpha values 
    vis = repmat(softSegments(:,:,1), [1 1 3]) .* repmat(colors(1,1,:), [size(softSegments, 1), size(softSegments, 2), 1]);
    for i = 2 : size(softSegments, 3)
        vis = vis + repmat(softSegments(:,:,i), [1 1 3]) .* repmat(colors(i,1,:), [size(softSegments, 1), size(softSegments, 2), 1]);    
    end
end

function colors = createPalette()
    % https://graphicdesign.stackexchange.com/questions/3682/where-can-i-find-a-large-palette-set-of-contrasting-colors-for-coloring-many-d
    % From P. Green-Armytage (2010): A Colour Alphabet and the Limits of Colour Coding. 
    % Colour: Design & Creativity (5) (2010): 10, 1-23
    % http://eleanormaclure.files.wordpress.com/2011/03/colour-coding.pdf
    colors = reshape(... 
                                    [0,117,220; 255,255,128; 43,206,72; 153,0,0; 128,128,128; 240,163,255; 153,63,0; 76,0,92;...
                                    0,92,49; 255,204,153; 148,255,181;...
                                    143,124,0; 157,204,0; 194,0,136; 0,51,128; 255,164,5;...
                                    255,168,187; 66,102,0; 255,0,16; 94,241,242; 0,153,143;...
                                    224,255,102; 116,10,255; 255,255,0; 255,80,5; 25,25,25]...
                                    , [26, 1, 3]) / 255;
    colors = [colors; colors];
end

function order = determineAlphaOrder(alphas)
    alphas = sum(alphas, 1);
    alphas = sum(alphas, 2);
    alphas = reshape(alphas, [size(alphas, 3) 1]);
    [~, order] = sort(alphas, 'descend');
end
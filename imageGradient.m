
% Yagiz Aksoy, 2016
% Implements H. Farid, E.P. Simoncelli, Differentiation of Discrete Multidimensional Signals, TIP 2004
% outColor switch makes the output 3-channel (color) or 1-channel, default true
% There are 3 different filtSize options (1, 4 or 6), which determines the number of taps in the filters
% and hence which neighborhood is used to approximate the derivatives

function [gradientMagnitude, gradientOrientation, xGradient, yGradient] = imageGradient(image, outColor, filtSize)
    % Set up variables legally
    if ~exist('outColor', 'var') || isempty(outColor)
        outColor = true;
    end
    if ~exist('filtSize', 'var') || isempty(filtSize)
        filtSize = 1;
    end
    if filtSize <= 3
        filtSize = 1;
    elseif filtSize <= 5
        filtSize = 4;
    else
        filtSize = 6;
    end
    image = im2double(image);

    % Set up one-dimensional filters
    switch filtSize
        case 1
            dk = [0.425287, -0.0000, -0.425287];
            kk = [0.229879, 0.540242, 0.229879];
        case 4
            dk = [0.0032, 0.0350, 0.1190, 0.1458, -0.0000, -0.1458, -0.1190, -0.0350, -0.0032];
            kk = [0.0009, 0.0151, 0.0890, 0.2349, 0.3201, 0.2349, 0.0890, 0.0151, 0.0009];
        case 6
            dk = [0.0001, 0.0019, 0.0142, 0.0509, 0.0963, 0.0878, 0.0000, -0.0878, -0.0963, -0.0509, -0.0142, -0.0019, -0.0001];
            kk = [0.0000, 0.0007, 0.0071, 0.0374, 0.1126, 0.2119, 0.2605, 0.2119, 0.1126, 0.0374, 0.0071, 0.0007, 0.0000];
    end

    % Repeat-pad image
    leftPad = image(:, 1, :);
    rightPad = image(:, end, :);
    image = [repmat(leftPad, [1 13 1]), image, repmat(rightPad, [1 13 1])];
    upPad = image(1, :, :);
    downPad = image(end, :, :);
    image  = [repmat(upPad, [13 1 1]); image; repmat(downPad, [13 1 1])];

    % Compute gradients
    yGradient = zeros(size(image));
    xGradient = zeros(size(image));
    for i = 1 : size(image, 3)
        yGradient(:,:,i) = conv2(dk, kk, image(:,:,i), 'same');
        xGradient(:,:,i) = conv2(kk, dk, image(:,:,i), 'same');
    end

    % Remove padding
    yGradient = yGradient(14 : end - 13, 14 : end - 13, :);
    xGradient = xGradient(14 : end - 13, 14 : end - 13, :);

    % Compute pixel-wise L2 norm if no color option is selected
    if ~outColor
        yGradient = sqrt(sum(yGradient .* yGradient, 3));
        xGradient = sqrt(sum(xGradient .* xGradient, 3));
    end

    % Compute polar-coordinate representation
    gradientMagnitude = sqrt(yGradient .* yGradient + xGradient .* xGradient);
    gradientOrientation = atan2(xGradient, yGradient);
end
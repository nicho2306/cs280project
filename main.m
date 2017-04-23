clear all, close all
% Input image
name = '8';
image = imread(['./pictures/' name '.jpg']);

% Process the image
processed = scanner(image, 1);

figure
imshow(image)
title('original image')
figure
imshow(processed);
title('scanned image')
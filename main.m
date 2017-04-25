clear all, close all
% Input image
name = '11';
image = imread(['./pictures/' name '.jpg']);

figure
imshow(image)
title('original image')

% Process the image
displayFigures = true;
processed = scanner(image, displayFigures);

figure
imshow(processed);
title('scanned image')
clear all, close all
% Input image
name = '11';
image = imread(['./pictures/' name '.jpg']);
%%
% Process the image
processed = scanner(image, 1);
processed = rgb2gray(processed)

figure
imshow(image)
title('original image')
figure
imshow(processed);
title('scanned image')
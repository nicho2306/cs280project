clear all, close all
% Input image
name = '09';
image = imread(['./pictures/' name '.jpg']);
%image = imresize(image, 0.5);
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
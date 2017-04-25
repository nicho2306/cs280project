clear all, close all
% Input image
name = '16';
image = imread(['./pictures/' name '.jpg']);
imshow(image);
title(name)
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
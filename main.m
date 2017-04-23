clear all, close all

name = '1';
image = imread(['./pictures/' name '.jpg']);

processed = scanner(image, 1);

figure
imshow(image)
title('original image')
figure
imshow(processed);
title('scanned image')
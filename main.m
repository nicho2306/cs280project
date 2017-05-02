clear all, close all
% input CSV file
filename = './pictures/cornersDataSet.csv';
corners_dataset = csvread(filename,1,1);
num_images = size(corners_dataset, 1) - 1;

%% Check corner detector accuracy over all dataset
good_counts = 0;
bad_images = [];
for i=0:num_images
    name = int2str(i);
    image = imread(['./pictures/' name '.jpg']);
    c_points = cornerDetector(image, false);
    d_points = corners_dataset(i+1,:);
    d_points = reshape(d_points, 2, 4);
    [m n p] = size(image);
    tol = (m + n)/2 * 0.02;
    good_pred = computeError(c_points, d_points, tol);
    if ~good_pred
        % store all images that don't get good corners
        bad_images = [bad_images i];
        % show the image that has bad corner prediction
        %displayCorners(image, c_points);
        %title(['image ', name, ' with bad corners'])

    end
    good_counts = good_counts + good_pred;
end

accuracy = good_counts/num_images

%% scanner: 
close all
i = 16; % choose an index for the image to load
name = int2str(i);
image = imread(['./pictures/' name '.jpg']);

% Process the image
displayFigures = true;
processed = scanner(image, displayFigures);

figure
imshow(processed);
title('scanned image')
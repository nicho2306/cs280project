clear all, close all
% input CSV file
filename = './pictures/cornersDataSet.csv';
corners_dataset = csvread(filename,1,1);
num_images = size(corners_dataset, 1) - 1;

%% Check corner detector accuracy over all dataset
good_counts = 0;
bad_images = [];
for i=18:num_images
    i
    name = int2str(i);
    image = imread(['./pictures/' name '.jpg']);
    c_points = cornerDetector(image, false);
    d_points = corners_dataset(i+1,:);
    d_points = reshape(d_points, 2, 4);
    good_pred = computeError(c_points, d_points, 20);
    if ~good_pred
        % store all images that don't get good corners
        bad_images = [bad_images i];
        % show the image that has bad corner prediction
        figure
        imshow(image)
        hold on
        plot(c_points(1,:),c_points(2,:),'rx','markersize', 30, 'linewidth',5)
        hold on
        plot(c_points(1,1:2),c_points(2,1:2),'Color','green','LineWidth',3)
        hold on
        plot(c_points(1,2:3),c_points(2,2:3),'Color','green','LineWidth',3)
        hold on
        plot(c_points(1,3:4),c_points(2,3:4),'Color','green','LineWidth',3)
        hold on
        plot([c_points(1,4),c_points(1,1)],[c_points(2,4),c_points(2,1)],...
            'Color','green','LineWidth',3)
        title(['image ', name, ' with bad corners'])
    end
    good_counts = good_counts + good_pred;
end

accuracy = good_counts/(num_images + 1)

%%
close all
i = 13; % choose an index for the image to load
name = int2str(i);
image = imread(['./pictures/' name '.jpg']);
cornerDetector(image, true);
%%
figure
imshow(image)
title('original image')

% Process the image
displayFigures = true;
processed = scanner(image, displayFigures);

figure
imshow(processed);
title('scanned image')
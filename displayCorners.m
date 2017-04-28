function [ ] = displayCorners(image, points )
% Display image and corner points
figure
imshow(image)
hold on
plot(points(1,:),points(2,:),'rx','markersize', 30, 'linewidth',5)
hold on
plot(points(1,1:2),points(2,1:2),'Color','green','LineWidth',3)
hold on
plot(points(1,2:3),points(2,2:3),'Color','green','LineWidth',3)
hold on
plot(points(1,3:4),points(2,3:4),'Color','green','LineWidth',3)
hold on
plot([points(1,4),points(1,1)],[points(2,4),points(2,1)],...
    'Color','green','LineWidth',3)
end


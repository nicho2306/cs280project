function target = scanner(image, displayFigures)
% Find corners in original image, and transform using homography
% The input arguments are:
%    image = original image to be transformed
%    angle = angle between camera and normal surface of the plane
% The output argument is:
%   target = transformed image


%% Find corners
points = cornerDetector(image, displayFigures);

%% Constraint the corner
lengthX = round(points(1,2) - points(1,3)) + 1;
lengthY = round(1.3*lengthX); %round(1.1*(points(2,2) - points(2,1))) + 1;
x4 = 1;
y4 = 1;
x1 = x4 + lengthX;
y1 = y4;
x3 = x4;
y3 = y4 + lengthY;
x2 = x1;
y2 = y3;
target_points = [x1 x2  x3 x4;
                 y1 y2  y3 y4];
                 
H = homography_solve(target_points, points);

%% Transform and crop the image
target = uint8(zeros(lengthY, lengthX, 3));
for x=1:size(target,2)
    for y=1:size(target,1)
        u = [x;y];
        v = homography_transform(u, H);
        target(y,x,:) = image(round(v(2)),round(v(1)),:);
    end
end
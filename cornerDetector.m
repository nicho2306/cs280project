function points = cornerDetector(image, displayCorners, displayHough)
I = rgb2gray(image);
I = imgaussfilt(I,2);
BW = edge(I, 'Canny', 0.5);

%% Hough transform
[H,T,R] = hough(BW);
P = houghpeaks(H,4,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(BW,T,R,P, 'FillGap',10000000000,'MinLength',1);

%% Find equation of lines
eq = zeros(3,length(lines));
for k=1:length(lines)
    X1 = lines(k).point1'; %(x1,y1)'
    X2 = lines(k).point2'; %(x2,y2)'
    eq(:,k) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
end

eq = zeros(3,length(lines));
sumX = zeros(1, 4);
sumY = zeros(1, 4);
for k=1:length(lines)
    X1 = lines(k).point1'; %(x1,y1)'
    X2 = lines(k).point2'; %(x2,y2)'
    sumX(k) = X1(1) + X2(1);
    sumY(k) = X1(2) + X2(2);
end
% The uppest line
i = find(sumY == min(sumY));
X1 = lines(i).point1'; %(x1,y1)'
X2 = lines(i).point2'; %(x2,y2)'
eq(:,1) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];

% The right line
i = find(sumX == max(sumX));
X1 = lines(i).point1'; %(x1,y1)'
X2 = lines(i).point2'; %(x2,y2)'
eq(:,2) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];

% The bottom line
i = find(sumY == max(sumY));
X1 = lines(i).point1'; %(x1,y1)'
X2 = lines(i).point2'; %(x2,y2)'
eq(:,3) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];

% The left line
i = find(sumX == min(sumX));
X1 = lines(i).point1'; %(x1,y1)'
X2 = lines(i).point2'; %(x2,y2)'
eq(:,4) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];

%% Find intersection of lines (corners)
points = zeros(2,4);   % first row: x-coordinate of the points, second row: y-coordinate
for i = 1:1:length(eq)
    
    if i ~= length(eq)
        l1 = eq(:, i);
        l2 = eq(:, i+1);
        intersect = cross(l1, l2);
        intersect = intersect/intersect(3);
        points(:,i) = intersect(1:2);
    else
        l1 = eq(:, i);
        l2 = eq(:, 1);
        intersect = cross(l1, l2);
        intersect = intersect/intersect(3);
        points(:,i) = intersect(1:2);
    end
end

%% Plot Hough lines & corners

if displayHough
    figure
    imshow(BW)
    hold on
    for k = 1:length(lines)
        xy = [lines(k).point1; lines(k).point2];
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
        % Plot beginnings and ends of lines
        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
    end
end

if displayCorners
    figure
    imshow(image)
    hold on
    plot(points(1,:),points(2,:),'x','markersize', 30, 'linewidth',5)
end
end
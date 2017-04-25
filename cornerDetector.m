function points = cornerDetector(image, displayCorners, displayHough)
% Find corners in original image
% The input arguments are:
%    image = original image to be filtered
%    displayCorners = whether to display image with plotted 4 corners
%    displayHough   = whether to display image with lines from Hough transform 
% The output argument is:
%    points = array of (x,y) coordinates of the 4 corners detected

%% convert to grayscale, lowpass filter, edge filter
I = rgb2gray(image);
I = imgaussfilt(I,2.5);
BW = edge(I, 'Canny', 0.3);

%% Hough transform
[H,T,R] = hough(BW);
P = houghpeaks(H,10,'threshold',ceil(0.3*max(H(:))));
LINES = houghlines(BW,T,R,P, 'FillGap',10000000000,'MinLength',1);

%% Find equation of lines
max_area = 0;
max_points = zeros(2,4);

% random sample 5000 times
for ii =1:5000
    % randomly get 4 candidate lines
    possibleK = 1:size(LINES,2);
    K = datasample(possibleK, 4,'Replace',false);
    lines = LINES(K);
    
    eq = zeros(3,length(lines));
    sumX = zeros(1, 4);
    sumY = zeros(1, 4);
    for k=1:4
        X1 = lines(k).point1'; %(x1,y1)'
        X2 = lines(k).point2'; %(x2,y2)'
        sumX(k) = X1(1) + X2(1);
        sumY(k) = X1(2) + X2(2);
    end
    
    % The top horizontal line
    i = find(sumY == min(sumY));
    i = i(randi(length(i)));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq(:,1) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    
    % The right vertical line
    i = find(sumX == max(sumX));
    i = i(randi(length(i)));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq(:,2) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    
    % The bottom horizontal line
    i = find(sumY == max(sumY));
    i = i(randi(length(i)));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq(:,3) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    
    % The left vertical line
    i = find(sumX == min(sumX));
    i = i(randi(length(i)));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq(:,4) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    
    
    % Find intersection of lines (corners) between the 4 random candidates
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
    
    % check whether these points are better candidates for corners
    % area cannot be too big, all points must be within the image pixels
    area = polyarea(points(1,:), points(2,:));
    if area > max_area && ~isnan(area) && ~any(any((points<0))) && ...
            ~any(points(1,:)>size(image,2)) && ~any(points(2,:)>size(image,1)) 
        max_area = area;
        max_points = points;
    end
    
end
points = max_points;

%% Plot Hough lines & corners

if displayHough
    figure
    imshow(BW)
    hold on
    for k = 1:length(LINES)
        xy = [LINES(k).point1; LINES(k).point2];
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
        % Plot beginnings and ends of lines
        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'markersize', 20, 'Color','yellow');
        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'markersize', 20, 'Color','red');
    end
end
%%
if displayCorners
    figure
    imshow(image)
    hold on
    plot(points(1,:),points(2,:),'rx','markersize', 30, 'linewidth',5)
   
end
end
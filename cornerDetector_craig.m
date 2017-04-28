function points = cornerDetector(image, displayFigures)
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
%I = medfilt2(I);
BW = edge(I, 'Canny', 0.3);

%% Hough transform
[H,T,R] = hough(BW);
P = houghpeaks(H,10,'threshold',ceil(0.3*max(H(:))));
LINES = houghlines(BW,T,R,P, 'FillGap',10000000000,'MinLength',1);

%% Find equation of lines


% random sample ~210 combinations of 10 lines
   
    eq = zeros(3, 4); % The outest lines
    eq2 = zeros(3, 4); % The second outest lines
    sumX = zeros(1, length(LINES));
    sumY = zeros(1, length(LINES));
    
    for k = 1:1:length(LINES)
        X1 = LINES(k).point1'; %(x1,y1)'
        X2 = LINES(k).point2'; %(x2,y2)'
        sumX(k) = X1(1) + X2(1);
        sumY(k) = X1(2) + X2(2);
    end
    
    lines = LINES;
    % The top horizontal line
    i = find(sumY == min(sumY));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq(:,1) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    % Second uppest line
    sumYY = sumY;
    sumYY(i) = [];
    i = find(sumYY == min(sumYY));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq2(:,1) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    
    % The right vertical line
    i = find(sumX == max(sumX));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq(:,2) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    % Second rightmost line
    sumXX = sumX;
    sumXX(i) = [];
    i = find(sumXX == max(sumXX));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq2(:,2) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    
    % The bottom horizontal line
    i = find(sumY == max(sumY));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq(:,3) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    % The second lowest line
    sumYY = sumY;
    sumYY(i) = [];
    i = find(sumYY == max(sumYY));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq2(:,3) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    
    % The left vertical line
    i = find(sumX == min(sumX));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq(:,4) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    % The second leftmost line
    sumXX = sumX;
    sumXX(i) = [];
    i = find(sumXX == min(sumXX));
    X1 = lines(i).point1'; %(x1,y1)'
    X2 = lines(i).point2'; %(x2,y2)'
    eq2(:,4) = [-(X2(2)-X1(2))/(X2(1)-X1(1)); 1; (X2(2)-X1(2))/(X2(1)-X1(1))*X1(1) - X1(2)];
    
    cell(1, :, :) = eq;
    cell(2, :, :) = eq2;
    allPoints = zeros(32, 4); % each two rows, and 4 columns is a combination of the points for the corners
    allAreas = zeros(1, 16); % 16 combinations in total
    
    % Find the combination of the two outest lines at the top, right,
    % bottom, and left. Calculate the area, and find the combination that
    % gives the biggest areas.
    num = 1;
    for a = 1:1:2
        for b = 1:1:2
            for c = 1:1:2
                for d = 1:1:2
                    Eq = [cell(a, :, 1)', cell(b, :, 2)', cell(c, :, 3)', cell(d, :, 4)'];
                    % Find intersection of lines (corners) between the 4 random candidates
                    points = zeros(2,4);   % first row: x-coordinate of the points, second row: y-coordinate
                    for i = 1:1:length(Eq)     
                        if i ~= length(Eq)
                        l1 = Eq(:, i);
                        l2 = Eq(:, i+1);
                        intersect = cross(l1, l2);
                        intersect = intersect/intersect(3);
                        points(:,i) = intersect(1:2);
                        else
                        l1 = Eq(:, i);
                        l2 = Eq(:, 1);
                        intersect = cross(l1, l2);
                        intersect = intersect/intersect(3);
                        points(:,i) = intersect(1:2);
                        end    
                    end
                    allPoints((num*2-1):(num*2),:) = points(:, :);
                    area = polyarea(points(1,:), points(2,:));
                    % Check if the points are out of bound
                    if  ~any(any((points<0))) && ~any(points(1,:)>size(image,2)) && ~any(points(2,:)>size(image,1)) 
                        allAreas(num) = area;
                    end
                    num = num + 1;
                end
            end
        end
    end
     i = find(allAreas == max(allAreas));
     points = allPoints((i*2-1):(i*2), :);
    


%% Plot Hough lines & corners

if displayFigures
    figure
    imshow(BW)
    title('edges of image')
    
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
    title('image with multiple line candidates')
    
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
    title('image with corners')
   
end
end
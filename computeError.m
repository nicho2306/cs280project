function good_prediction = computeError(cPoints, dPoints, r_threshold)
% Compute L-2 norm error between cornerDetector and data set points
% The input arguments are:
%    cPoints = (x,y) points from cornerDetector
%    dPoints = (x,y) points from dataset
%    r_threshold = margin radius of error per corner
% The output argument is:
%    good_prediction = boolean of whether predicted corners are good enough

distances = 999*ones(1,4);
for i=1:4
    distances(i) = sqrt((cPoints(1,i) - dPoints(1,i))^2 + (cPoints(2,i) - dPoints(2,i))^2);
end

good_prediction = ~any(distances > r_threshold);

end


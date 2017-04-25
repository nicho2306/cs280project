function error = computeError(cPoints, dPoints)
% Compute L-2 norm error between cornerDetector and data set points
% The input arguments are:
%    cPoints = (x,y) points from cornerDetector
%    dPoints = (x,y) points from dataset
% The output argument is:
%    error = L-2 norm error sqrt((x_c - x_d)^2 + (y_c - y_d)^2)

error = 0;

for i=1:4
    for j = 1:2
        error = error + ((cPoints(i,j) - dPoints(i,j))^2);
    end
end

error = sqrt(error);
end


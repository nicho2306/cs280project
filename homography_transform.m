% Author: Nicholas Sunjaya
% CS280
% HW1
% Problem 3.4

function v = homography_transform(u, H)
    U = [u;
        ones(1, size(u,2))];
    V = H*U;
    v = [V(1,:)./V(3,:);
         V(2,:)./V(3,:)];
end
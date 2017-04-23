% Author: Nicholas Sunjaya
% CS280
% HW1
% Problem 3.4

function H = homography_solve(u,v)
    N = size(v,2);
    V = reshape(v,[2*N,1]);
    A = [];
    for i=1:N
        A = [A;
             u(1,i) u(2,i)   1      0      0   0 -u(1,i)*v(1,i) -u(2,i)*v(1,i);
                  0      0   0 u(1,i) u(2,i)   1 -u(1,i)*v(2,i) -u(2,i)*v(2,i)];
    end
    h = A\V;
    H = [h(1) h(2) h(3);
         h(4) h(5) h(6);
         h(7) h(8)    1];
     
end
   
function x = fwdKin(q,xi)
% inputs: 
%       q is a 2x1 matrix where row(1)=theta 
%           and row(2)=phi joint space vectors: -pi<=theta<pi & 0<=phi<=pi
%       xi scalar value that determines positions along the arc (0-> base,
%           1-> tip)
%           eg. if you need the tip position for computing fwdKin, use xi=1.
%           For drawing arm, use xi values between 0 and 1.
% outputs:
%       x, 3x1 =[x y z] coordinates of the arc
x=zeros(3,1);
L=0.24;
theta = q(1);
phi = q(2);
if theta < -pi || theta >= pi
    msg = 'theta outside accepted bounds [-pi,pi). Given theta : ';
    error(strcat(msg, string(theta)))
end
if phi < 0 || phi > pi
    msg = 'phi outside accepted bounds [0, pi]. Given phi : ';
    error(strcat(msg, string(phi)))
end
x(1)=L * cos(theta) * xi ^ 2 * phi / 0.2e1 - L * cos(theta) * ... 
    xi ^ 4 * phi ^ 3 / 0.24e2 + L * cos(theta) * xi ^ 6 * phi ^ 5 / ...
    0.720e3 - L * cos(theta) * xi ^ 8 * phi ^ 7 / 0.40320e5 + L * cos(theta) * ...
    xi ^ 10 * phi ^ 9 / 0.3628800e7 - L * cos(theta) * xi ^ 12 * phi ^ 11 / 0.479001600e9 ...
    + L * cos(theta) * xi ^ 14 * phi ^ 13 / 0.87178291200e11;

x(2)=L * sin(theta) * xi ^ 2 * phi / 0.2e1 - L * sin(theta) * ...
    xi ^ 4 * phi ^ 3 / 0.24e2 + L * sin(theta) * xi ^ 6 * phi ^ 5 / ...
    0.720e3 - L * sin(theta) * xi ^ 8 * phi ^ 7 / 0.40320e5 + L * sin(theta) * ...
    xi ^ 10 * phi ^ 9 / 0.3628800e7 - L * sin(theta) * xi ^ 12 * phi ^ 11 / 0.479001600e9 ...
    + L * sin(theta) * xi ^ 14 * phi ^ 13 / 0.87178291200e11;

x(3)=L * xi - L * xi ^ 3 * phi ^ 2 / 0.6e1 + L * xi ^ 5 * phi ...
    ^ 4 / 0.120e3 - L * xi ^ 7 * phi ^ 6 / 0.5040e4 + L * xi ^ 9 * phi ^ 8 / ...
    0.362880e6 - L * xi ^ 11 * phi ^ 10 / 0.39916800e8 + L * xi ^ 13 * ...
    phi ^ 12 / 0.6227020800e10;
end

function q = invKin(X)
% inputs: 
%       X =[x;y] coordinates of the arc
% outputs:
%       q=[theta;phi] joint space vectors: -pi<=theta<pi & 0<=phi<=pi
L=0.24;
q=zeros(2,1);
x = X(1);
y = X(2);
if x == 0 && y ==0
    q(1) = 0;
    q(2) = 0;
elseif norm(X)>1e-6
    theta=atan2(y,x);
    phi=fzero(@(j) y - L / j * sin(theta) * (0.1e1 - cos(j)),pi/2);
    if isnan(theta) || isnan(phi)
        msg = 'fzero failed to resolve a solution. Given inputs [x y] : ';
        error(strcat(msg, mat2str([x y])))
    end
    q(1) = theta;
    q(2) = phi;
end
end

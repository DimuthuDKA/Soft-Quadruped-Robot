function C = gaitCalc(X, turnRadius, gaitType, legNumber)
%gaitCalc(gaitInput,turnRadius,gaitType, 1);

T = X(2);
x0 = X(3);
y0 = X(4);
r = X(5);
zeta = X(6);
t = mod(X(1),T);
L=0.24;

C = zeros(2,1);
%ground=r*cos(pi-zeta)+x0;

% Trot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmp(gaitType, 'Trotting')
    
    rRight=0.03;
    rLeft=0.05;
    
    %T2=(zeta/pi)*T;
    u=zeta/pi;
    %u=0.45;
    
    if legNumber == 2 || legNumber == 4
        t = t + (180/360)*T; % Flipping diagonal legs for the trotting phase.
    end
    t = mod(t,T);
    
    if legNumber == 1 || legNumber == 4
        if  t < u*T % 0.5*u*T is start location -> left edge of walk
            C(1) = x0+rRight*cos(zeta);
            C(2) = y0+rRight*sin(zeta)-(t/(u*T))*2*rRight*sin(zeta);
        elseif t >= u*T
            radians = (2*(pi-zeta)*(t-u*T))/(T-u*T);
            C(1) = rRight*sin(pi/2-radians-zeta)+x0;
            C(2) = -rRight*cos(pi/2-radians-zeta)+y0;
        end
    elseif legNumber == 2 || legNumber == 3
        if  t < u*T % 0.5*u*T is start location -> left edge of walk
            C(1) = x0+rLeft*cos(zeta);
            C(2) = y0+rLeft*sin(zeta)-(t/(u*T))*2*rLeft*sin(zeta);
        elseif t >= u*T
            radians = (2*(pi-zeta)*(t-u*T))/(T-u*T);
            C(1) = rLeft*sin(pi/2-radians-zeta)+x0;
            C(2) = -rLeft*cos(pi/2-radians-zeta)+y0;
        end
    elseif legNumber == 5
        C(1) = 0;
        C(2) = 0;
    end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TurnRight
if strcmp(gaitType, 'TurnRight')
    rRight=0.03;
    
    horizontalL = (2/3)*L;
    % adjust stride radius for right and left leg
    rRatio=(turnRadius + horizontalL)/(turnRadius - horizontalL);
    rLeft=rRatio*rRight;
    rLeft=0.05;
    u=zeta/pi;
    %u=0.5;
    
    if legNumber == 2 || legNumber == 4
        t = t + (180/360)*T; % Flipping diagonal legs for the trotting phase.
    end
    t = mod(t,T);
    
    if legNumber == 1 || legNumber == 4
        if  t < u*T % 0.5*u*T is start location -> left edge of walk
            C(1) = x0+rRight*cos(zeta);
            C(2) = y0+rRight*sin(zeta)-(t/(u*T))*2*rRight*sin(zeta);
        elseif t >= u*T
            radians = (2*(pi-zeta)*(t-u*T))/(T-u*T);
            C(1) = rRight*sin(pi/2-radians-zeta)+x0;
            C(2) = -rRight*cos(pi/2-radians-zeta)+y0;
        end
    elseif legNumber == 2 || legNumber == 3
        if  t < u*T % 0.5*u*T is start location -> left edge of walk
            C(1) = x0+rLeft*cos(zeta);
            C(2) = y0+rLeft*sin(zeta)-(t/(u*T))*2*rLeft*sin(zeta);
        elseif t >= u*T
            radians = (2*(pi-zeta)*(t-u*T))/(T-u*T);
            C(1) = rLeft*sin(pi/2-radians-zeta)+x0;
            C(2) = -rLeft*cos(pi/2-radians-zeta)+y0;
        end
        
    elseif legNumber == 5
        %phi = atan2(L,(horizontalL+turnRadius));
        % Adjust the curve of the body to match the turn angle.
        %psi = atan2(L,(horizontalL+turnRadius));
        psi = L/turnRadius;
        temp = fwdKin([-pi/2;psi],1);
        C(1) = temp(1);
        C(2) = temp(2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TurnLeft
if strcmp(gaitType, 'TurnLeft')
    rLeft=r;
    
    
    horizontalL = (2/3)*L;
    % adjust stride radius for right and left leg
    rRatio=(turnRadius + horizontalL)/(turnRadius - horizontalL);
    rRight=rRatio*rLeft;
    rRight=0.06;
    
    u=zeta/pi;
    u=0.5;
    
    if legNumber == 2 || legNumber == 4
        t = t + (180/360)*T; % Flipping diagonal legs for the trotting phase.
    end
    t = mod(t,T);
    
    if legNumber == 1 || legNumber == 4
        if  t < u*T % 0.5*u*T is start location -> left edge of walk
            C(1) = x0+rRight*cos(zeta);
            C(2) = y0+rRight*sin(zeta)-(t/(u*T))*2*rRight*sin(zeta);
        elseif t >= u*T
            radians = (2*(pi-zeta)*(t-u*T))/(T-u*T);
            C(1) = rRight*sin(pi/2-radians-zeta)+x0;
            C(2) = -rRight*cos(pi/2-radians-zeta)+y0;
        end
    elseif legNumber == 2 || legNumber == 3
        if  t < u*T % 0.5*u*T is start location -> left edge of walk
            C(1) = x0+rLeft*cos(zeta);
            C(2) = y0+rLeft*sin(zeta)-(t/(u*T))*2*rLeft*sin(zeta);
        elseif t >= u*T
            radians = (2*(pi-zeta)*(t-u*T))/(T-u*T);
            C(1) = rLeft*sin(pi/2-radians-zeta)+x0;
            C(2) = -rLeft*cos(pi/2-radians-zeta)+y0;
        end
        
    elseif legNumber == 5
        %phi = atan2(L,(horizontalL+turnRadius));
        % Adjust the curve of the body to match the turn angle.
        %psi = atan2(L,(horizontalL+turnRadius));
        psi = L/turnRadius;
        temp = fwdKin([pi/2;psi],1);
        C(1) = temp(1);
        C(2) = temp(2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Flip leg 1 and 2 due to opposing orientation
if legNumber == 2 || legNumber == 3
    C(2) = C(2) * -1;
end
if strcmp(gaitType, 'Clockwise') && (legNumber == 1 || legNumber == 4)
    C(1) = C(1) * -1;
elseif strcmp(gaitType, 'Anticlockwise') && (legNumber == 2 || legNumber == 3)
    C(1) = C(1) * -1;
end

C = C.';

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

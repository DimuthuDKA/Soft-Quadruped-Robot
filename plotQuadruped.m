%function fig = plotQuadruped(X,frontRightLegXYZ,frontLeftLegXYZ,bodyXYZ,backLeftLegXYZ,backRightLegXYZ,ground,gaitType,titleString)
function fig = plotQuadruped(X,ground,gaitType,titleString,viewOption)

animationResolution = 10; % Number of points to animate a single arm

animationFrontRightLegXYZ = zeros(3, animationResolution);
animationFrontLeftLegXYZ = zeros(3, animationResolution);
animationBackRightLegXYZ = zeros(3, animationResolution);
animationBackLeftLegXYZ = zeros(3, animationResolution);
animationBodyXYZ = zeros(3, animationResolution);

frontRightCurveParams = [X(1) X(2)].';
frontLeftCurveParams = [X(3) X(4)].';
backLeftCurveParams = [X(5) X(6)].';
backRightCurveParams = [X(7) X(8)].';
bodyCurveParams = [X(9) X(10)].';

psiVals = linspace(0,1,animationResolution);
for j = 1:animationResolution
    animationFrontRightLegXYZ(:,j) = fwdKin(frontRightCurveParams, psiVals(j));
    animationFrontLeftLegXYZ(:,j) = fwdKin(frontLeftCurveParams, psiVals(j));
    animationBodyXYZ(:,j) = fwdKin(bodyCurveParams, psiVals(j));
    animationBackLeftLegXYZ(:,j) = fwdKin(backLeftCurveParams, psiVals(j));
    animationBackRightLegXYZ(:,j) = fwdKin(backRightCurveParams, psiVals(j));
end

% Moving in X direction in robot coordinate frame.
% animationFrontRightLegXYZ = Rz(-pi/2)*animationFrontRightLegXYZ;
% animationFrontRightLegXYZ = Rx(pi/2)*animationFrontRightLegXYZ;
% 
% animationFrontLeftLegXYZ = Rz(pi/2)*animationFrontLeftLegXYZ;
% animationFrontLeftLegXYZ = Rx(-pi/2)*animationFrontLeftLegXYZ;
% 
% animationBodyXYZ = Rz(pi)*animationBodyXYZ;
% animationBodyXYZ = Ry(-pi/2)*animationBodyXYZ;
% tailOrigin = animationBodyXYZ(:,animationResolution);

animationBackRightLegXYZ = Rz(-pi/2)*animationBackRightLegXYZ;
animationBackRightLegXYZ = Rx(pi/2)*animationBackRightLegXYZ;

animationBackLeftLegXYZ = Rz(pi/2)*animationBackLeftLegXYZ;
animationBackLeftLegXYZ = Rx(-pi/2)*animationBackLeftLegXYZ;

%animationBodyXYZ = Rz(pi)*animationBodyXYZ;
animationBodyXYZ = Ry(pi/2)*animationBodyXYZ;
tailOrigin = animationBodyXYZ(:,animationResolution);

if strcmp(gaitType, "TurnRight")

animationFrontLeftLegXYZ = Rz(pi/2)*animationFrontLeftLegXYZ;
animationFrontLeftLegXYZ = Rx(-pi/2)*animationFrontLeftLegXYZ;
animationFrontLeftLegXYZ = Rz(-bodyCurveParams(2))*animationFrontLeftLegXYZ;
animationFrontLeftLegXYZ = animationFrontLeftLegXYZ + tailOrigin;

animationFrontRightLegXYZ = Rz(-pi/2)*animationFrontRightLegXYZ;
animationFrontRightLegXYZ = Rx(pi/2)*animationFrontRightLegXYZ;
animationFrontRightLegXYZ = Rz(-bodyCurveParams(2))*animationFrontRightLegXYZ;
animationFrontRightLegXYZ = animationFrontRightLegXYZ + tailOrigin;


% animationBodyXYZ = Rz(pi)*animationBodyXYZ;
% animationBodyXYZ = Ry(-pi/2)*animationBodyXYZ;
% bodyOrigin = animationBodyXYZ(:,100);


elseif strcmp(gaitType, "TurnLeft")
 
animationBackLeftLegXYZ = Rz(pi/2)*animationBackLeftLegXYZ;
animationBackLeftLegXYZ = Rx(-pi/2)*animationBackLeftLegXYZ;
animationBackLeftLegXYZ = Rz(-bodyCurveParams(2))*animationBackLeftLegXYZ;
animationBackLeftLegXYZ = animationBackLeftLegXYZ + tailOrigin;

animationBackRightLegXYZ = Rz(-pi/2)*animationBackRightLegXYZ;
animationBackRightLegXYZ = Rx(pi/2)*animationBackRightLegXYZ;
animationBackRightLegXYZ = Rz(-bodyCurveParams(2))*animationBackRightLegXYZ;
animationBackRightLegXYZ = animationBackRightLegXYZ + tailOrigin;

else
% animationFrontLeftLegXYZ = Rz(pi/2)*animationFrontLeftLegXYZ;
% animationBackLeftLegXYZ = Rx(-pi/2)*animationBackLeftLegXYZ;
% animationBackLeftLegXYZ = animationBackLeftLegXYZ + tailOrigin;
% 
% animationBackRightLegXYZ = Rz(-pi/2)*animationBackRightLegXYZ;
% animationBackRightLegXYZ = Rx(pi/2)*animationBackRightLegXYZ;
% animationBackRightLegXYZ = animationBackRightLegXYZ + tailOrigin;

animationFrontRightLegXYZ = Rz(-pi/2)*animationFrontRightLegXYZ;
animationFrontRightLegXYZ = Rx(pi/2)*animationFrontRightLegXYZ;
animationFrontRightLegXYZ = animationFrontRightLegXYZ + tailOrigin;

animationFrontLeftLegXYZ = Rz(pi/2)*animationFrontLeftLegXYZ;
animationFrontLeftLegXYZ = Rx(-pi/2)*animationFrontLeftLegXYZ;
animationFrontLeftLegXYZ = animationFrontLeftLegXYZ + tailOrigin;

end

% Moving in Y direction in robot coordinate frame.
% animationFrontRightLegXYZ = Ry(pi/2)*animationFrontRightLegXYZ;
% animationFrontLeftLegXYZ = Rz(pi)*animationFrontLeftLegXYZ;
% animationFrontLeftLegXYZ = Ry(-pi/2)*animationFrontLeftLegXYZ;
% animationBodyXYZ = Rz(-pi/2)*animationBodyXYZ;
% animationBodyXYZ = Rx(pi/2)*animationBodyXYZ;
% tailOrigin = animationBodyXYZ(:,animationResolution);
% animationBackLeftLegXYZ = Rz(pi)*animationBackLeftLegXYZ;
% animationBackLeftLegXYZ = Ry(-pi/2)*animationBackLeftLegXYZ;
% animationBackLeftLegXYZ = animationBackLeftLegXYZ + tailOrigin;
% animationBackRightLegXYZ = Ry(pi/2)*animationBackRightLegXYZ;
% animationBackRightLegXYZ = animationBackRightLegXYZ + tailOrigin;

fig = plot3(animationFrontRightLegXYZ(1,:), animationFrontRightLegXYZ(2,:),animationFrontRightLegXYZ(3,:),'b', ...
    animationFrontLeftLegXYZ(1,:), animationFrontLeftLegXYZ(2,:),animationFrontLeftLegXYZ(3,:),'m', ...
    animationBodyXYZ(1,:), animationBodyXYZ(2,:), animationBodyXYZ(3,:),'k', ...
    animationBackRightLegXYZ(1,:), animationBackRightLegXYZ(2,:), animationBackRightLegXYZ(3,:),'b',...
    animationBackLeftLegXYZ(1,:), animationBackLeftLegXYZ(2,:), animationBackLeftLegXYZ(3,:),'m'); 
hold on;
% 
% p1=plot3(frontRightLegXYZ(1,:),frontRightLegXYZ(2,:),frontRightLegXYZ(3,:)); 
% p2=plot3(frontLeftLegXYZ(1,:),frontLeftLegXYZ(2,:),frontLeftLegXYZ(3,:)); 
% p3=plot3(backLeftLegXYZ(1,:),backLeftLegXYZ(2,:),backLeftLegXYZ(3,:)); 
% p4=plot3(backRightLegXYZ(1,:),backRightLegXYZ(2,:),backRightLegXYZ(3,:)); 
% p1.Color = 'k'; p2.Color = 'k'; p3.Color = 'k'; p4.Color = 'k';
% p2.LineStyle = ':'; p2.LineWidth=1.5; p1.LineStyle = ':'; p1.LineWidth=1.5;
% p3.LineStyle = ':'; p3.LineWidth=1.5; p4.LineStyle = ':'; p4.LineWidth=1.5;

hold off;
fig(1).LineWidth = 8;
fig(2).LineWidth = 8;
fig(3).LineWidth = 8;
fig(4).LineWidth = 8;
fig(5).LineWidth = 8;

if viewOption == "Standard"
    % axis equal;
    title(titleString);
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    xlim([-0.1, 0.41]);
    ylim([-0.24, 0.24]);
    zlim([-ground, 0.05]);
    grid on;
    view(3);
elseif viewOption == "XY"
    title(titleString);
    xlabel("X");
    ylabel("Y");
    xlim([-0.1, 0.41]);
    ylim([-0.24, 0.24]);
    zlim([-ground, 0.05]);
    grid on;
    view([0 0 1]);
elseif viewOption == "XZ"
    title(titleString);
    xlabel("X");
    zlabel("Z");
    xlim([-0.1, 0.41]);
    ylim([-0.24, 0.24]);
    zlim([-ground, 0.05]);
    grid on;
    view([0 -1 0]);
elseif viewOption == "YZ"
    title(titleString);
    xlabel("Y");
    zlabel("Z");
    xlim([-0.1, 0.41]);
    ylim([-0.24, 0.24]);
    zlim([-ground, 0.05]);
    grid on;
    view(90,0);
end

hold off;
fig = gcf;
end
clear;
generateAnimation = 1;
makeGraphs = 0;
writeCurveParameterFiles = 0;
plotsingleanimation=0;
close all;

%rho= (L1/pphi)*(1-cos(pphi)); %in meters

% pphi= (pi/180)*39.75;
L1=0.24;
% rho=0.06;
% xd=rho+(1-cos(pphi))*L1/pphi;
% d= (L1/pphi)*sin(pphi);

animationResolution = 10; % Number of points to animate a single arm.
n=100; % Number of points to calculate for one period.
T = n; % Period length


r = 0.03; % Stride radius
x0 = 0; %Origind
turnRadius= 0.5;

y0 = -0.00000001; %Origin
zeta = 0;
%zeta = 2*pi/2-0.00000000000001;

% ground=r*cos(pi-zeta)+x0;
ground=r*cos(zeta)+x0;
del = 0.08; % Set the time between frames in the generated gifs. Units in seconds

frontRightCurveParams = zeros(2,n);
frontLeftCurveParams = zeros(2,n);
backLeftCurveParams = zeros(2,n);
backRightCurveParams = zeros(2,n);
bodyCurveParams = zeros(2,n);

frontRightLegXY = zeros(2,n);
frontLeftLegXY = zeros(2,n);
backLeftLegXY = zeros(2,n);
backRightLegXY = zeros(2,n);
bodyXY = zeros(2,n);
D= zeros (1,n);

originalFrontRightLegXY = zeros(2,n);
originalFrontLeftLegXY = zeros(2,n);
originalBackLeftLegXY = zeros(2,n);
originalBackRightLegXY = zeros(2,n);
originalBodyXY = zeros(2,n);

times = linspace(1,n);
% for gait = [Trotting, TurnLeft, TurnRight]
for gaitType = "TurnRight" % Change this for the desired gait
    
    if ~exist(strcat(pwd,"\generatedGifs\",string(gaitType),"_Gait\"), 'dir')
        mkdir(strcat(pwd,"\generatedGifs\",string(gaitType),"_Gait\"))
    end
    if ~exist(strcat(pwd,"\generatedCSVs\",string(gaitType),"_Gait\"), 'dir')
        mkdir(strcat(pwd,"\generatedCSVs\",string(gaitType),"_Gait\"))
    end
    if ~exist(strcat(pwd,"\generatedGraphs\",string(gaitType),"_Gait\"), 'dir')
        mkdir(strcat(pwd,"\generatedGraphs\",string(gaitType),"_Gait\"))
    end
    
    dataFileMatrix = zeros(n, 11);
    
    for i = 1:n
        gaitInput = [i T x0 y0 r zeta];
        
        originalFrontRightLegXY(:,i) = gaitCalc(gaitInput,turnRadius,gaitType, 1);
        originalFrontLeftLegXY(:,i) = gaitCalc(gaitInput,turnRadius,gaitType, 2);
        originalBackLeftLegXY(:,i) = gaitCalc(gaitInput,turnRadius,gaitType, 3);
        originalBackRightLegXY(:,i) = gaitCalc(gaitInput,turnRadius,gaitType, 4);
        originalBodyXY(:,i) = gaitCalc(gaitInput,turnRadius,gaitType, 5);
        
        frontRightCurveParams(:,i) = invKin(originalFrontRightLegXY(:,i));
        frontLeftCurveParams(:,i) = invKin(originalFrontLeftLegXY(:,i));
        backLeftCurveParams(:,i) = invKin(originalBackLeftLegXY(:,i));
        backRightCurveParams(:,i) = invKin(originalBackRightLegXY(:,i));
        bodyCurveParams(:,i) = invKin(originalBodyXY(:,i)); %zeros(2,1);
        
        dataFileMatrix(i,:) = [i,frontRightCurveParams(:,i).',frontLeftCurveParams(:,i).',backLeftCurveParams(:,i).',backRightCurveParams(:,i).',bodyCurveParams(:,i).'];
        
                temp = fwdKin(frontRightCurveParams(:,i),1);
                frontRightLegXY(:,i) = temp(1:2,:);
                temp = fwdKin(frontLeftCurveParams(:,i),1);
                frontLeftLegXY(:,i) = temp(1:2,:);
                temp = fwdKin(backLeftCurveParams(:,i),1);
                backLeftLegXY(:,i) = temp(1:2,:);
                temp = fwdKin(backRightCurveParams(:,i),1);
                backRightLegXY(:,i) = temp(1:2,:);
                temp = fwdKin(zeros(2,1),1);
                bodyXY(:,i) = temp(1:2,:);
        
%         temp = fwdKin(frontRightCurveParams(:,i),1);
%         frontRightLegXYZ(:,i) = temp(1:3,:);
%         temp = fwdKin(frontLeftCurveParams(:,i),1);
%         frontLeftLegXYZ(:,i) = temp(1:3,:);
%         temp = fwdKin(backLeftCurveParams(:,i),1);
%         backLeftLegXYZ(:,i) = temp(1:3,:);
%         temp = fwdKin(backRightCurveParams(:,i),1);
%         backRightLegXYZ(:,i) = temp(1:3,:);
%         temp = fwdKin(bodyCurveParams(:,i),1);
%         bodyXYZ(:,i) = temp(1:3,:);
        
        
        % Generate Animation
        if generateAnimation == 1
            curveParams = [frontRightCurveParams(:,i) frontLeftCurveParams(:,i) backLeftCurveParams(:,i) backRightCurveParams(:,i) bodyCurveParams(:,i)];
            myFigure = plotQuadruped(curveParams,ground,gaitType,strcat(string(gaitType)," - Standard View"), "Standard");
            %plotQuadruped(X,ground,titleString,viewOption)
            %plotQuadruped(X,frontRightLegXYZ,frontLeftLegXYZ,bodyXYZ,backLeftLegXYZ,backRightLegXYZ,ground,gaitType,titleString)
            
            frame = getframe(myFigure);
            
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if i == 1
                imwrite(imind,cm,strcat(pwd,"\generatedGifs\",string(gaitType),"_Gait\",string(gaitType),'_StandardView','_Animation.gif'),'gif', 'Loopcount',inf,'DelayTime',del);
            else
                imwrite(imind,cm,strcat(pwd,"\generatedGifs\",string(gaitType),"_Gait\",string(gaitType),'_StandardView','_Animation.gif'),'gif','WriteMode','append','DelayTime',del);
            end
            clf('reset')
            
            myFigure = plotQuadruped(curveParams,ground,gaitType,strcat(string(gaitType)," - XY View"), "XY");
            frame = getframe(myFigure);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if i == 1
                imwrite(imind,cm,strcat(pwd,"\generatedGifs\",string(gaitType),"_Gait\",string(gaitType),'_XY_View','_Animation.gif'),'gif', 'Loopcount',inf,'DelayTime',del);
            else
                imwrite(imind,cm,strcat(pwd,"\generatedGifs\",string(gaitType),"_Gait\",string(gaitType),'_XY_View','_Animation.gif'),'gif','WriteMode','append','DelayTime',del);
            end
            clf('reset')
            
            myFigure = plotQuadruped(curveParams,ground,gaitType,strcat(string(gaitType)," - XZ View"),"XZ");
            frame = getframe(myFigure);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if i == 1
                imwrite(imind,cm,strcat(pwd,"\generatedGifs\",string(gaitType),"_Gait\",string(gaitType),'_XZ_View','_Animation.gif'),'gif', 'Loopcount',inf,'DelayTime',del);
            else
                imwrite(imind,cm,strcat(pwd,"\generatedGifs\",string(gaitType),"_Gait\",string(gaitType),'_XZ_View','_Animation.gif'),'gif','WriteMode','append','DelayTime',del);
            end
            clf('reset')
            
            myFigure = plotQuadruped(curveParams,ground,gaitType,strcat(string(gaitType)," - YZ View"), "YZ");
            frame = getframe(myFigure);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File
            if i == 1
                imwrite(imind,cm,strcat(pwd,"\generatedGifs\",string(gaitType),"_Gait\",string(gaitType),'_YZ_View','_Animation.gif'),'gif', 'Loopcount',inf,'DelayTime',del);
            else
                imwrite(imind,cm,strcat(pwd,"\generatedGifs\",string(gaitType),"_Gait\",string(gaitType),'_YZ_View','_Animation.gif'),'gif','WriteMode','append','DelayTime',del);
            end
            clf('reset')
            
        end
    end
    
    % The curve parameters data file is a csv file organized by
    %(time, FRTheta, FRPhi, FLTheta, FLPhi, BLTheta, BLPhi, BRTheta, BRPhi, BodyTheta, BodyPhi)
    if writeCurveParameterFiles == 1
        writematrix(dataFileMatrix, strcat(pwd,"\generatedCSVs\",string(gaitType),"_Gait\",string(gaitType),"_curveParams.csv"))
    end
    
    % All these graphs are saved into separate image files
    if makeGraphs == 1
        close all force
        n=100;
        f=0.25;T=1/f;N=3;r=0.018;
        t=linspace(0,T,n+1)';t(end)=[];
        times=t;
        
        r1 = [0.4660 0.6740 0.1880]; %green
        r2 = [0.6350 0.0780 0.1840]; %blood red
        r3 = 'r';%[0.8500 0.3250 0.0980];
        r4 = [0 0.4470 0.7410];%blue
        
        
        figure(1);
        %subplot(2,1,1);
        plot(times, frontRightLegXYZ(1,:), 'r','LineWidth',2); hold on;
        plot(times, frontLeftLegXYZ(1,:), 'b','LineWidth',2);
        plot(times, backLeftLegXYZ(1,:), 'g','LineWidth',2);
        plot(times, backRightLegXYZ(1,:), 'm','LineWidth',2);
        grid;
        fig1 = gca;
        xlabel("Time [s]");
        ylabel("X [m]");
        set(gca,'FontSize',24); %set(gca,'FontWeight','bold');
        set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
        %         xlim([0 4]);
        %         xticks([0:1:4]);
        %         ylim([0.03 0.15]);
        %         yticks([0.03:0.03:0.15]);
        hold off;
        
        
        figure(2);
        %subplot(2,1,2);
        plot(times, frontRightLegXYZ(2,:), 'r','LineWidth',2); hold on;
        plot(times, frontLeftLegXYZ(2,:), 'b','LineWidth',2);
        plot(times, backLeftLegXYZ(2,:), 'g','LineWidth',2);
        plot(times, backRightLegXYZ(2,:), 'm','LineWidth',2);
        grid;
        fig2 = gca;
        set(gca,'FontSize',24); %set(gca,'FontWeight','bold');
        set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
        %         lgd = legend('Leg1','Leg2','Leg3','Leg4');%lgd.NumColumns = 3;
        %         lgd.FontSize = 8;
        xlabel("Time [s]");
        ylabel("Y [m]");
        %         xlim([0 4]);
        %         xticks([0:1:4]);
        %         yticks([-0.06:0.03:0.06]);
        %         ylim([-0.06 0.06]);
        hold off;
        
        
        figure(3);
        plot(times, frontRightLegXYZ(3,:),'Color',r1,'LineWidth',3); hold on;
        plot(times, frontLeftLegXYZ(3,:),'Color',r2,'LineWidth',3);
        plot(times, backLeftLegXYZ(3,:),'Color',r3,'LineWidth',3);
        plot(times, backRightLegXYZ(3,:),'Color',r4,'LineWidth',3);
        grid;
        fig3 = gca;
        set(gca,'FontSize',22); %set(gca,'FontWeight','bold');
        set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
        %         lgd = legend('Leg1','Leg2','Leg3','Leg4');%lgd.NumColumns = 3;
        %         lgd.FontSize = 8;
        xlabel("Time [s]");
        ylabel("Z [m]");
                 xlim([0 4]);
                 xticks([0:1:4]);
                 yticks([0.12:0.02:0.26]);
                 ylim([0.12 0.26]);
        hold off;
        
        figure(4);
        plot(t, frontRightCurveParams(1,:), 'r','LineWidth',2); hold on;
        plot(t, frontLeftCurveParams(1,:), 'b','LineWidth',2);
        plot(t, backLeftCurveParams(1,:), 'g','LineWidth',2);
        plot(t, backRightCurveParams(1,:), 'm','LineWidth',2);
        grid;
        %title("Curve Parameters");
        xlabel("Time [s]");
        %ylabel("Theta [rad]");
        fig2 = gca;
        set(gca,'FontSize',24); %set(gca,'FontWeight','bold');
        set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
        %         xlim([0 4]);
        %         xticks([0:1:4]);
        %         yticks([-0.6:0.2:0.6]);
        %         ylim([-0.6 0.6]);
        hold off;
        
        
        figure(5);
        plot(t, frontRightCurveParams(2,:), 'r','LineWidth',2); hold on;
        plot(t, frontLeftCurveParams(2,:), 'b','LineWidth',2);
        plot(t, backLeftCurveParams(2,:), 'g','LineWidth',2);
        plot(t, backRightCurveParams(2,:), 'm','LineWidth',2);
        grid;
        xlabel("Time [s]");
        %ylabel("Phi [rad]");
        set(gca,'FontSize',24); %set(gca,'FontWeight','bold');
        set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
        %         xlim([0 4]);
        %         xticks([0:1:4]);
        %         yticks([0:0.4:2.2]);
        %         ylim([0 2.2]);
        hold off;
        
        figure(6);
        subplot(2,1,1);
        plot(originalFrontRightLegXY(2,:), -originalFrontRightLegXY(1,:), 'blue');
        grid;
        title("Original Gait Shape - Front Right");
        xlabel("Y");
        ylabel("X");
        axis equal
        
        temp = zeros(3, n);
        for i = 1:n
            temp(:,i) = fwdKin(frontRightCurveParams(:,i), 1);
        end
        subplot(2,1,2);
        plot(temp(2,:), -temp(1,:), 'red');
        grid;
        title("Forward Kinematic Gait Shape - Front Right");
        xlabel("Y");
        ylabel("X");
        axis equal;
        fig3 = gca;
        
        figure(7);
        subplot(2,1,1);
        plot(originalFrontLeftLegXY(2,:), -originalFrontLeftLegXY(1,:), 'blue');
        grid;
        title("Original Gait Shape - Front Left");
        xlabel("Y");
        ylabel("X");
        axis equal
        
        temp = zeros(3, n);
        for i = 1:n
            temp(:,i) = fwdKin(frontLeftCurveParams(:,i), 1);
        end
        subplot(2,1,2);
        plot(temp(2,:), -temp(1,:), 'red');
        grid;
        title("Forward Kinematic Gait Shape - Front Left");
        xlabel("Y");
        ylabel("X");
        axis equal;
        fig3 = gca;
        
        figure(8);
        i=1;
        curveParams = [frontRightCurveParams(:,i) frontLeftCurveParams(:,i) backLeftCurveParams(:,i) backRightCurveParams(:,i) bodyCurveParams(:,i)];
        myFigure = plotQuadruped(curveParams,frontRightLegXYZ,frontLeftLegXYZ,bodyXYZ,backLeftLegXYZ,backRightLegXYZ,ground,gaitType,strcat(string(gaitType)));
        frame = getframe(myFigure); grid on; axis equal;
        %plot3(frontRightLegXYZ(1,:),frontRightLegXYZ(2,:),frontRightLegXYZ(3,:)); grid on;
        fig8 = gca;
        set(gca,'FontSize',16); %set(gca,'FontWeight','bold');
        set(gca,'xcolor','k');set(gca,'ycolor','k');set(gca,'zcolor','k');
        
        %for fig = [fig1,fig2,fig3]
        %            saveas(fig,strcat(pwd,"\generatedGraphs\",string(gaitType),"_Gait\",fig.Title.String,".png"));
        %end
        
    end
end
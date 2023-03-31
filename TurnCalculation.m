    L=0.24;
    turnRadius=0.6;
r=0.04;
rRight=r;
    
    horizontalL = (2/3)*L;
    % adjust stride radius for right and left leg
    rRatio=(turnRadius + 0.12)/(turnRadius - 0.16);
    rLeft=rRatio*rRight
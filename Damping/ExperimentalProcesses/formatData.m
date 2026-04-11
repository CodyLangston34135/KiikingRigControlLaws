function [theta,thetaDot,blueDist,redDist,time] = formatData(theta,blueDist,redDist,time,datum,maxChange,debugTheta)
    
    indexN = 1:length(time);
    if debugTheta
        figure(1)
        subplot(3,1,1)
        plot(indexN,theta)
    end
    
    %% Remove Dropped Points
    % Dropped points for red distance
    [~, index] = find(redDist~=0);
    theta = theta(index);
    blueDist = blueDist(index);
    redDist = redDist(index);
    time = time(index);

    % Dropped points for blue distance
    [~, index] = find(blueDist~=0);
    theta = theta(index);
    blueDist = blueDist(index);
    redDist = redDist(index);
    time = time(index);

    % Inaccurate points for red distance
    [~, index] = find(and(redDist>=(datum(1)-maxChange(2)),redDist<=(datum(2)+maxChange(2))));
    theta = theta(index);
    blueDist = blueDist(index);
    redDist = redDist(index);
    time = time(index);

    % Innacurate points for blue distance
    [~, index] = find(and(blueDist>=(datum(3)-maxChange(3)),blueDist<=(datum(3)+maxChange(3))));
    theta = theta(index);
    blueDist = blueDist(index);
    redDist = redDist(index);
    time = time(index);

    % Dropped points for theta
    [~, index] = find(theta~=0);
    theta = theta(index);
    blueDist = blueDist(index);
    redDist = redDist(index);
    time = time(index);
    thetaDot = (theta(2:end)-theta(1:end-1))./(time(2:end)-time(1:end-1));

    if debugTheta
        figure(1)
        subplot(3,1,2)
        plot(index,theta)

        figure(2)
        subplot(2,1,1)
        plot(index(2:end),thetaDot)
    end

    [~, index] = find(and(thetaDot<=maxChange(1),thetaDot>=-maxChange(1)));
    theta = theta(index);
    blueDist = blueDist(index);
    redDist = redDist(index);
    time = time(index);
    thetaDot = (theta(2:end)-theta(1:end-1))./(time(2:end)-time(1:end-1));

    thetaDotDot = (thetaDot(2:end)-thetaDot(1:end-1))./(time(2:end-1)-time(1:end-2));
    [~, index] = find(and(thetaDotDot<=10000,thetaDotDot>=-10000));
    theta = theta(index);
    blueDist = blueDist(index);
    redDist = redDist(index);
    time = time(index);
    thetaDot = (theta(2:end)-theta(1:end-1))./(time(2:end)-time(1:end-1));

    if debugTheta
        figure(1)
        subplot(3,1,3)
        plot(index,theta)

        figure(2)
        subplot(2,1,2)
        plot(index(2:end),thetaDot)

        figure(3)
        subplot(2,1,1)
        plot(index,redDist)
        subplot(2,1,2)
        plot(index,blueDist)
    end
    
    thetaDot = gradient(theta,time);
end
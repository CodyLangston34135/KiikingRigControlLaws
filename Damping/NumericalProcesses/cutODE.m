function [theta,thetaDot,overflow] = cutODE(theta,thetaDot,endCond,numCycle,odeTime)
    
    % Store Variables
    if endCond(2)==3
        tFinal = moveTime;
    else
        tFinal = odeTime(1);
    end
    tOverflow = odeTime(2);
    
    % Check when ODE initially goes over condition
    if endCond(2)==1        % Checks for theta
        thetaF = endCond(1);
        tempCond = or(and(theta(2:end)>=thetaF,theta(1:end-1)<=thetaF),and(theta(2:end)<=thetaF,theta(1:end-1)>=thetaF));
        tempCond = find(tempCond);
        for i = 1:numCycle
            [endPoint, ind] = min(tempCond);
            tempCond(ind) = [];
        end
    elseif endCond(2)==2    % Checks for thetaDot
        thetaDotF = endCond(1);
        tempCond = or(and(thetaDot(2:end)>=thetaDotF,thetaDot(1:end-1)<=thetaDotF),and(thetaDot(2:end)<=thetaDotF,thetaDot(1:end-1)>=thetaDotF));
        tempCond = find(tempCond);
        endPoint = min(tempCond);
        for i = 1:numCycle
            [endPoint, ind] = min(tempCond);
            tempCond(ind) = [];
        end
    else
        endPoint = length(theta)-1;
    end

    % Check if ODE did not meet condition
    overflow = 0;
    if tFinal >= tOverflow          % Overflow, will never meet condition
        disp('ODE overflow')
        if abs(theta(end))<=pi
            theta(end) = 0;
            thetaDot(end) = 0;
            overflow = 1;           % ODE overflow where system goes to 0
        else
            overflow = 2;           % ODE overflow where system goes over
        end
    elseif size(endPoint,1)==0      % If ODE did not meet condition, run recursively
        [thetaTemp,thetaDotTemp,time,overflow] = runODE(theta(end),thetaDot(end),endCond,mu,g,swing,state,moveTime,[tFinal*2 tOverflow],options);
        theta(end+1:end+length(thetaTemp)) = thetaTemp;
        thetaDot(end+1:end+length(thetaTemp)) = thetaDotTemp;
    else
        theta = theta(1:endPoint+1);
        thetaDot = thetaDot(1:endPoint+1);

        if endCond(2)==1
            thetaDotF = thetaDot(end-1)+(endCond(1)-theta(end-1))*(thetaDot(end)-thetaDot(end-1))/(theta(end)-theta(end-1));
            theta(end) = endCond(1);
            thetaDot(end) = thetaDotF;
        elseif endCond(2)==2
            thetaF = theta(end-1)+(endCond(1)-thetaDot(end-1))*(theta(end)-theta(end-1))/(thetaDot(end)-thetaDot(end-1));
            theta(end) = thetaF;
            thetaDot(end) = endCond(1);
        end
    end
end
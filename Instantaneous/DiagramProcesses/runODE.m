function [theta,thetaDot,time,overflow] = runODE(thetaI,thetaDotI,endCond,numCycle,mu,g,swing,state,moveTime,odeTime,options)

    %% Store Variables
    if endCond(2)==3
        tFinal = moveTime;
    else
        tFinal = odeTime(1);
    end
    tOverflow = odeTime(2);

    % Initial Conditions
    IC = [thetaI, thetaDotI];
    tSpan = [0 tFinal];

    % Properties
    Ms = swing(1);
    Ls = swing(2);
    Hs = swing(3);
    Is = swing(4);
    Mp = state(1,1);

    % Inital Position
    LpI = state(1,2);
    HpI = state(1,3);
    IpI = state(1,4);
    
    % Final Position
    LpF = state(2,2);
    HpF = state(2,3);
    IpF = state(2,4);

%% Calculate ODE Parameters

    % Time-Varying Properties
    Ldot = (LpF-LpI)/moveTime;
    Lp = @(t) LpI + Ldot*t;
    Ip = @(t) Mp*Lp(t)^2;

    % Coefficients
    I = @(t) Is+Ip(t);
    ML = @(t) Mp*Lp(t) + Ms*Ls/2;
    Cf = @(t) 2*Mp*Lp(t)*Ldot;
        
    % ODE Setup
    thDot = @(t,th) [th(2); (-Cf(t)/I(t)*(th(2))-mu*sign(th(2))-g*ML(t)/I(t)*sin(th(1)))];
    
    %% Run ODE
    [time, temp] = ode45(thDot,tSpan,IC,options);
    
    theta = temp(:,1);
    thetaDot = temp(:,2);

    %% Cut ODE to Correct Size
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
        [thetaTemp,thetaDotTemp,time,overflow] = runODE(theta(end),thetaDot(end),endCond,numCycle,mu,g,swing,state,moveTime,[tFinal*2 tOverflow],options);
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
 function [theta,thetaDot,time] = runODE(thetaI,thetaDotI,moveInput,swing,state,odeTime,endCond)

    % Initial Conditions
    IC = [thetaI, thetaDotI];
    tSpan = [0 odeTime];

    % Properties
    % Global
    mu = moveInput(1);
    g = moveInput(2);
    moveTime = moveInput(3);

    % Swing
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
    Mnc = @(t) Mp*Lp(t);
    Cf = @(t) 2*Mp*Lp(t)*Ldot;
        
    % ODE Setup
    % Conservative
    % thDot = @(t,th) [th(2); (-Cf(t)/I(t)*(th(2))-mu*sign(th(2))-g*ML(t)/I(t)*sin(th(1)))];
    thDot = @(t,th) [th(2); (-Cf(t)/I(t)*(th(2))-mu*th(2)-g*ML(t)/I(t)*sin(th(1)))];

    % ODE Options
    switch endCond(2)
        case 3
            options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3);
        case 4
            options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3,'Events',@(time,thTemp)myControlEvent(time,thTemp,endCond));
        case 5
            options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3,'Events',@(time,thTemp)myControlEvent(time,thTemp,endCond));
        case 6
            global storeThTemp storeTime
            storeThTemp = [thetaI,thetaDotI];
            storeTime = 1;
            options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3,'Events',@(time,thTemp)myAccelEvent(time,thTemp,endCond));
        case 7
            options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3,'Events',@(time,thTemp)myExperimentalControlEvent(time,thTemp,endCond));
        case 8
            options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3,'Events',@(time,thTemp)myExperimentalControlEvent(time,thTemp,endCond));
        otherwise
            options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3,'Events',@(time,thTemp)myEvent(time,thTemp,endCond));
    end
    
    %% Run ODE
    [time, thTemp] = ode45(thDot,tSpan,IC,options);
    
    theta = thTemp(:,1);
    thetaDot = thTemp(:,2);
end

function [value, isterminal, direction] = myEvent(time, thTemp, endCond)
    value      = [thTemp(endCond(2))-endCond(1), time-2];
    isterminal = [1 1];   % Stop the integration
    direction  = [0 0];
end

function [value, isterminal, direction] = myControlEvent(time, thTemp, endCond)
    value      = [thTemp(endCond(2)-3)-endCond(1), time-2];
    isterminal = [1 1];   % Stop the integration
    direction  = [0 0];
end

function [value, isterminal, direction] = myAccelEvent(time, thTemp, endCond)
    global storeThTemp storeTime
    dThdt = abs((thTemp(1)-storeThTemp(1))/(time(1)-storeTime(1)));
    value      = [dThdt-endCond(1), time-2];
    isterminal = [1 1];   % Stop the integration
    direction  = [0 0];
    storeThTemp = thTemp;
    storeTime = time;
end

function [value, isterminal, direction] = myExperimentalControlEvent(time, thTemp, endCond)
    value      = [thTemp(endCond(2)-6)-endCond(1), time-2];
    isterminal = [1 1];   % Stop the integration
    direction  = [endCond(3) 0];
end

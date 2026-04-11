function [theta,thetaDot,tau] = runNDODE(thetaI,thetaDotI,state,moveInput,endCond)

% Properties
sigma = state(1);
delta = state(2);
beta = state(3);
zeta = state(4);

% Initial Conditions
IC = [thetaI, thetaDotI];
tSpan = [0 beta];

switch moveInput
    case 1       
        omegaNtau = @(tau) sqrt((3+6*sigma*(1-delta/beta*tau))/(3+6*sigma*(1-delta/beta*tau)^2)*g);
        coefDisp = 1;
        coefVel = @(tau) -6*sigma*delta/(beta*(1+3*sigma*(1-delta/beta*tau))) + 2*zeta;
        thDot = @(t,th) [th(2); -coefVel(t)*th(2)-coefDisp*sin(th(1))];
    case 2
        omegaNtau = @(tau) sqrt((3+6*sigma*(1-delta*(1-tau/beta)))/(2+6*sigma*(1-delta*(1-tau/beta))^2));
        coefDisp = 1;
        coefVel = @(tau) 6*sigma*delta/(beta*(1+3*sigma*(1-delta*(1-tau/beta)))) + 2*zeta;
        thDot = @(t,th) [th(2); -coefVel(t)*th(2)-coefDisp*sin(th(1))];
end

switch endCond(2)
    case 1
        options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3,'Events',@(time,thTemp)myEvent(time,thTemp,endCond));
    case 2
        options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3,'Events',@(time,thTemp)myEvent(time,thTemp,endCond));
    case 3
        options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3);
end

[tauTemp, thTemp] = ode45(thDot,tSpan,IC,options);
 
tau = tauTemp;
theta = thTemp(:,1);
thetaDot = thTemp(:,2);

end

function [value, isterminal, direction] = myEvent(time, thTemp, endCond)
    value      = [thTemp(endCond(2))-endCond(1), time-2];
    isterminal = [1 1];   % Stop the integration
    direction  = [0 0];
end
%% ODE Test
clc; clear all; close all;

%% Nondimensional Parameters
% Tau = 2*pi/omegaN
sigma = 2.5;            % Nondimensional Mass Mr/Ms
delta = 0.05;            % Nondimensional Length 1-Lst/Ls
beta = 01*2*pi;               % Nondimensional Time Tst/Tperiod*2*pi
zeta = 0.12;            % Nondimensional Damping

thetaI = pi/8;
thetaDotI = pi;

% Energy
% nrgKEsq = @(thetaDot,sigma) (2+3*sigma)/(2+6*sigma)*thetaDot^2;
% nrgPEsq = @(theta) 1-cos(theta);
% nrgTsq = @(theta,thetaDot,sigma) nrgKEsq(thetaDot,sigma)+nrgPEsq(theta); 

nrgTsq = @(theta,thetaDot) 1/6*thetaDot^2+1/3*(1-cos(theta));

%% Physical Parameters
Ms = 0.15;
Ls = 0.25;
g = 9.81;

omegaNt = sqrt((3+6*sigma)/(2+6*sigma)*g/Ls);
omegaN = sqrt(g*((2*sigma+1)/2*Ms*Ls)/((1+3*sigma)/3*Ms*Ls^2));
Lst = (1-delta)*Ls;
Lsq = Ls;
Mp = sigma*Ms;
moveTime = beta/omegaN;

Is = 1/3*Ms*Ls^2;   % Moment of inertia swing
Isq = Mp*Lsq^2;     % Moment of inertia squat
Ist = Mp*Lst^2;     % Moment of inertia stand

datumPE = Ls; % PE Datum
Hs = datumPE-Ls/2;
Hsq = datumPE-Lsq;
Hst = datumPE-Lst;
Swing = [Ms Ls Hs Is];
Squat = [Mp Lsq Hsq Isq];
Stand = [Mp Lst Hst Ist];

% Energy
nrgKEr = @(thetaDot,I) 1/2*I*thetaDot.^2;               % Rotational Kinetic Energy
nrgPE = @(theta,M,L,H) M*g*(H+L*(1-cos(theta)));        % Gravitational Potential Energy
nrgT = @(theta,thetaDot,Ip,Lp,Hp) (nrgPE(theta,Mp,Lp,Hp)+nrgPE(theta,Ms,Ls/2,Hs))...    % Potential
     + (nrgKEr(thetaDot,Ip)+nrgKEr(thetaDot,Is));                                       % Rotational Kinetic


%% Nondimensional ODE Standing
% Time-Varying Parameters
% omega2 = @(tau) sqrt((3+6*sigma*(1-delta*tau/beta))/(2+6*sigma*(1-delta*tau/beta)^2));
% coefDisp = @(tau) omega2(tau)^2;
% coefVel = @(tau) 2*zeta*omega2(tau)+(6*sigma*delta*(delta*tau/beta-1))/(beta*(1+3*sigma*(1-delta*tau/beta)^2));

omegaNtau = @(tau) sqrt((3+6*sigma*(1-delta/beta*tau))/(3+6*sigma*(1-delta/beta*tau)^2)*g);
coefDisp = 1;
coefVel = @(tau) -6*sigma*delta/(beta*(1+3*sigma*(1-delta/beta*tau))) + 2*zeta;

thDot = @(t,th) [th(2); -coefVel(t)*th(2)-coefDisp*sin(th(1))];

% ODE settings
endCond = [0 1];
IC = [thetaI, 0];
tSpan = [0 beta];
options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3,'Events',@(time,thTemp)myEvent(time,thTemp,endCond));
[tauTemp, thTemp] = ode45(thDot,tSpan,IC);

% Results
tauNDSt = tauTemp;
timeNDSt = tauTemp/omegaN;
thetaNDSt = thTemp(:,1);
thetaDotNDSt = thTemp(:,2)*omegaN;

nrgIND = nrgTsq(thetaI,0)
nrgFND = nrgTsq(thTemp(end,1),thTemp(end,2))
%% Physical ODE Standing
% Time-Varying Properties
Ldot = (Lst-Lsq)/moveTime;
Lp = @(t) Lsq + Ldot*t;
Ip = @(t) Mp*Lp(t)^2;
I = @(t) Is+Ip(t);
ML = @(t) Mp*Lp(t) + Ms*Ls/2;
omegaNp = @(t) sqrt(g*ML(t)/I(t));
Mnc = @(t) Mp*Lp(t);
Cf = @(t) 2*Mp*Lp(t)*Ldot;
thDot = @(t,th) [th(2); (-Cf(t)/I(t)*(th(2))-2*zeta*omegaNp(t)*th(2)-g*ML(t)/I(t)*sin(th(1)))];

% ODE settings
endCond = [0 1];
IC = [thetaI, 0];
tSpan = [0 beta/omegaN];
[timeSt, thTemp] = ode45(thDot,tSpan,IC);

% Results
thetaSt = thTemp(:,1);
thetaDotSt = thTemp(:,2);

nrgI = nrgT(thetaI,0,Isq,Lsq,Hsq);
nrgF = nrgT(thetaSt(end),thetaDotSt(end),Isq,Lsq,Hsq);

%% Standing Plots
figure(1)
subplot(2,2,1)
plot(thetaSt,thetaDotSt);
hold on
plot(thetaNDSt,thetaDotNDSt,'--');
xlabel('Theta')
ylabel('ThetaDot')
title('Standing')

figure(1)
subplot(2,2,2)
plot(timeSt*omegaN,thetaSt)
hold on
plot(timeNDSt*omegaN,thetaNDSt,'--')
xlabel('Time')
ylabel('Theta')
title('Standing')

figure(2)
subplot(1,2,1)
plot(tauNDSt,thetaNDSt)
title('Standing')

%% Nondimensional ODE Squatting
% Time-Varying Parameters
% omega4 = @(tau) sqrt((3+6*sigma*(1-delta*(1-tau/beta)))/(2+6*sigma*(1-delta*(1-tau/beta))^2));
% coefDisp = @(tau) omega4(tau)^2;
% coefVel = @(tau) 2*zeta*omega4(tau)+(6*sigma*(1-delta^2*(1-tau/beta)))/(beta*(1+3*sigma*(1-delta*(1-tau/beta))^2));

omegaNtau = @(tau) sqrt((3+6*sigma*(1-delta*(1-tau/beta)))/(2+6*sigma*(1-delta*(1-tau/beta))^2));
coefDisp = 1;
coefVel = @(tau) 6*sigma*delta/(beta*(1+3*sigma*(1-delta*(1-tau/beta)))) + 2*zeta;

thDot = @(t,th) [th(2); -coefVel(t)*th(2)-coefDisp*sin(th(1))];

% ODE settings
endCond = [0 2];
IC = [0, thetaDotI/omegaN];
tSpan = [0 beta];
options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3,'Events',@(time,thTemp)myEvent(time,thTemp,endCond));
[tauTemp, thTemp] = ode45(thDot,tSpan,IC);

% Results
tauNDSq = tauTemp;
timeNDSq = tauTemp/omegaN;
thetaNDSq = thTemp(:,1);
thetaDotNDSq = thTemp(:,2)*omegaN;

%% Physical ODE Squatting
% Time-Varying Properties
Ldot = (Lsq-Lst)/moveTime;
Lp = @(t) Lst + Ldot*t;
Ip = @(t) Mp*Lp(t)^2;
I = @(t) Is+Ip(t);
ML = @(t) Mp*Lp(t) + Ms*Ls/2;
Mnc = @(t) Mp*Lp(t);
Cf = @(t) 2*Mp*Lp(t)*Ldot;
thDot = @(t,th) [th(2); (-Cf(t)/I(t)*(th(2))-2*zeta*omegaN*th(2)-g*ML(t)/I(t)*sin(th(1)))];

% ODE settings
endCond = [0 2];
IC = [0, thetaDotI];
tSpan = [0 beta/omegaN];
[timeSq, thTemp] = ode45(thDot,tSpan,IC);

% Results
thetaSq = thTemp(:,1);
thetaDotSq = thTemp(:,2);

% moveInput = [mu g moveTime];
% [tempTheta,tempThetaDot,~] = runODE(0,pi,moveInput,Swing,[Stand; Squat],moveTime,[0 2]);

%% Squatting Plots
figure(1)
subplot(2,2,3)
plot(thetaSq,thetaDotSq);
hold on
plot(thetaNDSq,thetaDotNDSq,'--');
xlabel('Theta')
ylabel('ThetaDot')
title('Squatting')

figure(1)
subplot(2,2,4)
plot(timeSq,thetaSq)
hold on
plot(timeNDSq,thetaNDSq,'--')
xlabel('Time')
ylabel('Theta')
title('Squatting')

figure(2)
subplot(1,2,2)
plot(tauNDSq,thetaNDSq)
title('Squatting')

function [value, isterminal, direction] = myEvent(time, thTemp, endCond)
    value      = [thTemp(endCond(2))-endCond(1), time-2];
    isterminal = [1 1];   % Stop the integration
    direction  = [0 0];
end



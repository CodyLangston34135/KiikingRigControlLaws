%% ODE45Testing
% Cody Langston
clc; clear all; close all;

L = 0.3;
g = 9.8;
m = 0.5;
H = 0;

mu = 0;

Ms = 1;
Ls = 0.3;

I = m*L^2;
Is = 1/3*Ms*Ls^2;

nrgKEr = @(thetaDot,I) 1/2*I*thetaDot.^2;               % Rotational Kinetic Energy
nrgKEt = @(thetaDot,M,L) 1/2*M*(thetaDot*L).^2;         % Translational Kinetic Energy
nrgPE = @(theta,M,L,H) M*g*(H+L*(1-cos(theta)));        % Gravitational Potential Energy


A = I+Is;
B = m*L + Ms*Ls/2;

thDot = @(t,th) [th(2); -B*g*t*sin(th(1))/A];
% thDot = @(t,th) [th(2); -g/L*sin(th(1))];
tSpan = [0 3];
IC = [pi/6,0];
options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3);
[time, thetaOut] = ode45(thDot,tSpan, IC, options);
plot(thetaOut(:,1),thetaOut(:,2));



PE1 = nrgPE(pi/6,m,L,H)+nrgPE(pi/6,Ms,Ls/2,Ls/2);
PE2 = nrgPE(0,m,L,H)+nrgPE(0,Ms,Ls/2,Ls/2);

sqrt((PE1-PE2)/(1/2*(Is+I)))


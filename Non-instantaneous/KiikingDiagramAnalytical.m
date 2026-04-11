%% Kiiking Analytical Diagrams
% Cody Langston
clc; clear all; close all;
addpath(fullfile(pwd, 'NumericalProcesses'));
% pool = parpool("Threads");

%% User Inputs
% Parameters
Ms = 0.157;     % Mass of swing
Mp = 0.3;       % Mass of person
Ls = 0.305;     % Length of swing
Lsq = 0.305;    % Length of squat
Lst = 0.28;   % Length of stand
g = 9.8;
Damping = .01;


% Moment of inertia
Is = 1/3*Ms*Ls^2;   % Moment of inertia swing
Isq = Mp*Lsq^2;     % Moment of inertia squat
Ist = Mp*Lst^2;     % Moment of inertia stand

% Combinations
MLsq = (Ms*Ls/2+Mp*Lsq);
MLst = (Ms*Ls/2+Mp*Lst);
MIsq = Is+Isq;
MIst = Is+Ist;

% Compile States
Swing = [Ms Ls Is];
Squat = [Mp Lsq Isq];
Stand = [Mp Lst Ist];

omegaN = sqrt(g*(Ms*Ls/2+Mp*Lsq)/(Is+Isq));
mu = 2*Damping*omegaN;

% Initialize Variables
nStep = 100;
n = 150;
theta = zeros(0,1);
thetaDot = zeros(0,1);

%% Energy Calculations

nrgKE = @(thetaDot,I) 1/2*I*thetaDot.^2;               % Rotational Kinetic Energy
nrgPE = @(theta,M,L) M*g*(L*(1-cos(theta)));        % Gravitational Potential Energy
nrgDampPE = @(theta,x,Mp,mu) -mu*sqrt(2*g*MLsq/MIsq)*trapz(x,(cos(x)-cos(theta)).^(1/2));
nrgDampKE = @(thetaDot,x,Mp,mu) -mu*trapz(x,sqrt(1/2*MIst*thetaDot^2-g*MLst*(1-cos(x))));
% Total energy
nrgT = @(theta,thetaDot,Ip,Lp) (nrgPE(theta,Mp,Lp)+nrgPE(theta,Ms,Ls/2))...    % Potential
     + (nrgKE(thetaDot,Ip)+nrgKE(thetaDot,Is));                                       % Rotational Kinetic

% Create Analytical Expressions
calcThetaDot2 = @(theta1,Mp,mu,x) sqrt((nrgPE(theta1,Ms,Ls/2)+nrgPE(theta1,Mp,Lsq)-abs(nrgDampPE(theta1,x,Mp,mu)))/(1/2*MIsq));
calcThetaDot3 = @(thetaDot2,Mp,mu) thetaDot2*MIsq/MIst;
calcNRG4 = @(thetaDot3,thetaEnd,Mp,mu,x) (nrgKE(thetaDot3,Is)+nrgKE(thetaDot3,Isq)-abs(nrgDampKE(thetaDot3,x,Mp,mu)));
calcNRGEnd = @(theta1,Mp,mu) calcNRG4(calcThetaDot3(calcThetaDot2(theta1)),theta1);

%% Calculate State
theta1 = pi;
x = linspace(0,theta1,100);
nrgOver = nrgPE(theta1,Ms,Ls/2)+nrgPE(theta1,Mp,Lst);
thetaDot2 = calcThetaDot2(theta1,Mp,mu,x);
thetaDot3 = calcThetaDot3(thetaDot2,Mp,mu);
nrg4 = calcNRG4(thetaDot3,pi,Mp,mu,x);

if nrg4>=nrgOver && thetaDot2 >= 0 && thetaDot3 >= 0
    % disp('Goes over')
    state = pi;
else
    error = 1;
    theta = [pi 0];
    while error > 0.1
        thetaTemp = (theta(1)+theta(2))/2;
        nrgStart = nrgPE(thetaTemp,Ms,Ls/2)+nrgPE(thetaTemp,Mp,Lst);

        x = linspace(0,thetaTemp,100);
        thetaDot2Temp = calcThetaDot2(thetaTemp,Mp,mu,x);
        thetaDot3Temp = calcThetaDot3(thetaDot2Temp,Mp,mu);
        nrgTemp = calcNRG4(thetaDot3Temp,thetaTemp,Mp,mu,x);

        if nrgTemp>nrgStart && thetaDot2Temp>=0 && thetaDot3Temp >=0
            theta(2) = thetaTemp;
        else
            theta(1) = thetaTemp;
        end

        error = (theta(1)-theta(2))/pi;
    end

    if thetaTemp <= 0.2
        % disp('Damps Out')
        state = 0;
    else
        % disp('Limit Cycle Oscillation')
        state = thetaTemp;
    end
end

% surf(Mass,Damping,storeState)
% shading interp
% 
% xlabel('Mass');
% ylabel('Damping');
% zlabel('Theta')




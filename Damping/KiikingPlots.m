%% Kiiking Diagrams
% Cody Langston
clc; clear all; close all;
addpath(fullfile(pwd, 'NumericalProcesses'));
% pool = parpool("Threads");

%% User Inputs
% Nondimensional Inputs
sigma = 4;            % Nondimensional Mass Mr/Ms
delta = 0.05;            % Nondimensional Length 1-Lst/Ls
beta = 0.1*2*pi;        % Nondimensional Time Tst/Tperiod*2*pi
zeta = 0.02;             % Nondimensional Damping
thetaI = linspace(0,5*pi/6,100);

% zeta = linspace(0,0.15,20);
% beta = linspace(0.025,0.3,20)*2*pi;
delta = linspace(0,0.2,20);

% Energy Calculations
nrgTsq = @(theta,thetaDot) 1/6*thetaDot^2+1/3*(1-cos(theta));

for j = 1:size(delta,2)
    state = [sigma, delta(j), beta, zeta];
    parfor i = 1:size(thetaI,2)
        [thetaStore, thetaDotStore, tauStore, stateStore, thetaOptim(i)] = runCycle(thetaI(i),state);
    
        nrgI = nrgTsq(thetaI(i),0);
        nrgF = nrgTsq(thetaStore(end),0);
        nrgDel(i) = nrgF-nrgI;
        thetaDel(i) = abs(thetaStore(end))-thetaI(i);
    end

    nrgStore(j,:) = nrgDel;
    thetaStore(j,:) = thetaDel;
    thetaOptimStore(j,:) = thetaOptim;

end

surf(thetaI,delta,nrgStore)
% surf(thetaI,beta/(2*pi),nrgStore)
colormap jet
colorbar
clim([-0.05 0.05]);
shading interp
xlabel('Initial Theta (rad)')

figure(2)
contour(thetaI,delta,nrgStore,[-0.01 0.01],"ShowText",true)
% contour(thetaI,beta/(2*pi),nrgStore,[-0.01 0 0.01],"ShowText",true)
xlabel('Initial Theta (rad)')

figure(3)
surf(thetaI,delta,thetaOptimStore)
colormap jet
colorbar
shading interp
xlabel('Initial Theta (rad)')

figure(4)
contour(thetaI,delta,thetaOptimStore,"ShowText",true)
xlabel('Initial Theta (rad)')
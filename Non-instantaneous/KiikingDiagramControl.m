%% Kiiking Control Diagrams
% Cody Langston
clc; clear all; close all;
addpath(fullfile(pwd, 'NumericalProcesses'));
addpath(fullfile(pwd, 'Figures'));

%% User Inputs
% Experiment
experimentalFile = 'output12.mat';

% Parameters
Ms = 0.157;     % Mass of swing
Mp = 0.3;       % Mass of person
Ls = 0.406;     % Length of swing
Lsq = 0.305;    % Length of squat
Lst = 0.29;     % Length of stand
g = 9.8;        % Gravity
Damping = .0165;     % Damping Coefficent

% Controls
velN = 0.1;          % Velocity threshold
accelN = 0;        % Acceleration threshold
dispN = 0.001;               % Displacement threshold

datumPE = Ls; % PE Datum
Hs = datumPE-Ls/2;
Hsq = datumPE-Lsq;
Hst = datumPE-Lst;

% Moment of inertia
Is = 1/3*Ms*Ls^2;   % Moment of inertia swing
Isq = Mp*Lsq^2;     % Moment of inertia squat
Ist = Mp*Lst^2;     % Moment of inertia stand

% Initial Conditions
thetaI = -0.6;              % Initial angle (rad)
thetaDotI = -3.3;              % Initial velocity (rad/s)

% Time parameters
odeTime = 2;          % Simulation and overflow time for ODE45
nStep = 150;                % Number of total cycles
n = 50;                    % Number of undamped values per cycle
moveTime = 0.1;
numCycle = 1;

% Compile States
Swing = [Ms Ls Hs Is];
Squat = [Mp Lsq Hsq Isq];
Stand = [Mp Lst Hst Ist];

omegaN = sqrt(g*(Ms*Ls/2+Mp*Lsq)/(Is+Isq));
mu = 2*Damping*omegaN;

% Initialize Variables
thetaStep = thetaI;
thetaDotStep = thetaDotI;
nrgStep = [0 0 0];                                      % [PEsq=KEsq=PEi KEst=PEst PEf]
endCond = [0 0];                                        % [endValue endType] 1 = theta, 2 = thetaDot, 3 = moveTime, 4 = displacement, 5 = velocity, 6 = accleration
overflow = 0;                                           % 0 = no overflow, 1 = damps out, 2 = goes over
stage = 0;
storeState = zeros(1,0);
storeTime = zeros(1,0);
storeTheta = zeros(1,0);
storeThetaDot = zeros(1,0);
currentTime = 0;

%% Energy Calculations
nrgKEr = @(thetaDot,I) 1/2*I*thetaDot.^2;               % Rotational Kinetic Energy
nrgPE = @(theta,M,L,H) M*g*(H+L*(1-cos(theta)));        % Gravitational Potential Energy

% Total energy
nrgT = @(theta,thetaDot,Ip,Lp,Hp) (nrgPE(theta,Mp,Lp,Hp)+nrgPE(theta,Ms,Ls/2,Hs))...    % Potential
     + (nrgKEr(thetaDot,Ip)+nrgKEr(thetaDot,Is));                                       % Rotational Kinetic

%% Initial Conditions
if thetaDotI==0 || thetaDotI*thetaI<0
    nrgTotal(1) = nrgT(thetaI,thetaDotI,Isq,Lsq,Hsq);
    stage = 1;
elseif thetaI ==0 || thetaDotI*thetaI>0
    nrgTotal(1) = nrgT(thetaI,thetaDotI,Ist,Lst,Hst);
    stage = 2;
end

tic
for i = 1:nStep

    %% Gaining Kineteic Energy
    if stage == 1 % Gaining kinetic energy need to squat

        % Ending System Properties
        endCond = [dispN 4];
        % endCond = [accelN 6];
        moveInput = [mu g moveTime];
        [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Squat; Squat],odeTime,endCond);   % Stationary Damped

        [thetaMove,thetaDotMove,timeMove] = runODE(thetaDamp(end),thetaDotDamp(end),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Moving Damped

        % Adjust time
        timeDamp = timeDamp+currentTime;
        currentTime = timeDamp(end);
        timeMove = timeMove+currentTime;
        currentTime = timeMove(end);

        % Phase Diagram
        figure(1)
        subplot(1,2,1)
        hold on
        plot(thetaDamp,thetaDotDamp,'g');
        plot(thetaMove,thetaDotMove,'r');

        % Theta Diagram
        figure(2)
        subplot(2,2,1)
        hold on
        plot(timeDamp,thetaDamp,'g')
        plot(timeMove,thetaMove,'r')
        subplot(2,2,3)
        hold on
        plot(timeDamp,thetaDotDamp,'g')
        plot(timeMove,thetaDotMove,'r')
        
        % Store Values
        thetaStep(i+1,:) = real(thetaMove(end));
        thetaDotStep(i+1,:) = real(thetaDotMove(end));
        storeState(end+1:end+size(timeDamp,1)-1,1) = 0; % Stationary
        storeState(end+1:end+size(timeMove,1)-1,1) = 2; % Squatting
        storeTime(end+1:end+size(timeDamp,1)-1,1) = timeDamp(2:end);
        storeTime(end+1:end+size(timeMove,1)-1,1) = timeMove(2:end);
        storeTheta(end+1:end+size(thetaDamp,1)-1,1) = thetaDamp(2:end);
        storeTheta(end+1:end+size(thetaMove,1)-1,1) = thetaMove(2:end);
        storeThetaDot(end+1:end+size(thetaDotDamp,1)-1,1) = thetaDotDamp(2:end);
        storeThetaDot(end+1:end+size(thetaDotMove,1)-1,1) = thetaDotMove(2:end);

        nrgTotal(i+1) = nrgT(thetaMove(end),thetaDotMove(end),Ist,Lst,Hst);
        nrgCheck = nrgTotal(i+1) - nrgT(0,0,Ist,Lst,Hst);
        stage = 2;

    %% Gaining Potential Energy
    elseif stage == 2 % Gaining potential energy need to stand
        
        % Ending System Properties
        endCond = [velN,5];
        moveInput = [mu g moveTime];
        [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Stand; Stand],odeTime,endCond);   % Stationary Damped

        [thetaMove,thetaDotMove,timeMove] = runODE(thetaDamp(end),thetaDotDamp(end),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Moving Damped

        % Adjust time
        timeDamp = timeDamp+currentTime;
        currentTime = timeDamp(end);
        timeMove = timeMove+currentTime;
        currentTime = timeMove(end);

        % Phase Diagram
        figure(1)
        subplot(1,2,1)
        hold on
        plot(thetaDamp,thetaDotDamp,'g');
        plot(thetaMove,thetaDotMove,'b');

        % Theta Diagram
        figure(2)
        subplot(2,2,1)
        hold on
        plot(timeDamp,thetaDamp,'g')
        plot(timeMove,thetaMove,'b')
        subplot(2,2,3)
        hold on
        plot(timeDamp,thetaDotDamp,'g')
        plot(timeMove,thetaDotMove,'b')

        % Store Values
        thetaStep(i+1,:) = real(thetaMove(end));
        thetaDotStep(i+1,:) = real(thetaDotMove(end));
        storeState(end+1:end+size(timeDamp,1)-1,1) = 0; % Stationary
        storeState(end+1:end+size(timeMove,1)-1,1) = 1; % Standing
        storeTime(end+1:end+size(timeDamp,1)-1,1) = timeDamp(2:end);
        storeTime(end+1:end+size(timeMove,1)-1,1) = timeMove(2:end);
        storeTheta(end+1:end+size(thetaDamp,1)-1,1) = thetaDamp(2:end);
        storeTheta(end+1:end+size(thetaMove,1)-1,1) = thetaMove(2:end);
        storeThetaDot(end+1:end+size(thetaDotDamp,1)-1,1) = thetaDotDamp(2:end);
        storeThetaDot(end+1:end+size(thetaDotMove,1)-1,1) = thetaDotMove(2:end);

        nrgTotal(i+1) = nrgT(thetaMove(end),thetaDotMove(end),Isq,Lsq,Hsq);
        nrgCheck = nrgTotal(i+1)-nrgT(0,0,Isq,Lsq,Hsq);
        stage = 1;
    end

    

    if abs(thetaStep(i+1,1))>=pi
        disp('Goes Over')
        overflow = 2;

        figure(1)
        subplot(1,2,1)
        xline(pi,'--k')
        xline(-pi,'--k')
        xlim([-pi-0.5 pi+0.5])
        break
    elseif nrgCheck <= 0.002
        disp('Damps Out')
        break
    elseif i == 75
        disp('Limit Cycle Oscillation')
        break
    end
end

%% Format Plots
figure(1)
subplot(1,2,1)
title('Phase Diagram')
xlabel('Angular Displacement (rad)')
ylabel('Angular Velocity (rad/s)')

figure(2)
subplot(2,2,1)
title('Angular Displacement')
xlabel('Time (s)')
ylabel('Angular Displacement (rad)')
subplot(2,2,3)
title('Angular Velocity')
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')

%% Create Energy Plots
figure
hold on
plot(nrgTotal(2:2:end),'b')
plot(nrgTotal(1:2:end),'r')
legend('Max Potential Energy','Max Kinetic Energy');
xlabel('Half-Oscillations')

nrgOver = nrgT(pi,0,Ist,Lst,Hst);
if overflow == 2
    yline(nrgOver,'--k');
    legend('Energy Gained','Energy Lost','Required Energy')
else
    legend('Energy Gained','Energy Lost')
end

%% Create Time Lag Embed Diagram

timeN = 0:0.005:currentTime;
thetaN = interp1(storeTime,storeTheta,timeN);
thetaDotN = interp1(storeTime,storeThetaDot,timeN);
stateN = interp1(storeTime,storeState,timeN);
stateN = ceil(stateN);

for i = 2:length(stateN)-1
    if stateN(i-1) == 0 && stateN(i) == 1 && stateN(i+1) == 2
        stateN(i) = 2;
    elseif stateN(i-1) == 2 && stateN(i) == 1 && stateN(i+1) == 0
        stateN(i) = 2;
    end
end

% figure
% plot(thetaN)
periodSteps = 30;

% Plot time lag embed plot
figure(4)
subplot(1,2,1)
for i = 2:size(thetaN,2)-periodSteps-1
    hold on
    if stateN(i)==1
        plot([thetaN(i), thetaN(i+1)],[thetaN(i+periodSteps), thetaN(i+periodSteps+1)],'b');
    elseif stateN(i)==2
        plot([thetaN(i), thetaN(i+1)],[thetaN(i+periodSteps), thetaN(i+periodSteps+1)],'r');
    else
        plot([thetaN(i), thetaN(i+1)],[thetaN(i+periodSteps), thetaN(i+periodSteps+1)],'g');
    end
end
title('Time Lag Embed Diagram')
ylabel('Angular Displacement (rad)')
xlabel('Angular Displacement (rad)')
xlim([-4 4])
ylim([-4 4])

plotExperimental(experimentalFile)

toc
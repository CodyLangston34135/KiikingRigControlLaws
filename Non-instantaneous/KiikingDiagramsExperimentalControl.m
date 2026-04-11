%% Kiiking Experimental Control Diagrams
clc; clear all; close all;

addpath(fullfile(pwd, 'NumericalProcesses'));
addpath(fullfile(pwd, 'ExperimentalProcesses'));
addpath(fullfile(pwd, 'ExperimentalData'));
addpath(fullfile(pwd, 'Figures'));

dataInput = 'output14.mat';
load(dataInput) % timeN thetaN thetaDotSN thetaDotDotSN distN moveIndex

%% Format Standing/Squatting times
index = 1;
storeRun = [1,1,moveIndex(1)];
for i = 1:size(moveIndex,2)-1
    if moveIndex(i) == moveIndex(i+1)
        storeRun(index,2) = i+1;
    else
        index=index+1;
        storeRun(index,:) = [i+1, i+1, moveIndex(i+1)];
    end
end

for i = 1:size(storeRun,1)
    if storeRun(i,3) == 0
        distBounds(i) = mean(distN(storeRun(i,1):storeRun(i,2)));
    end
end
moveTime = (storeRun(:,2)-storeRun(:,1))*1/60;

%% User Inputs
% Parameters
Ms = 0.157;     % Mass of swing
Mp = 0.3;       % Mass of person
Ls = 0.3;     % Length of swing
Lp = distBounds*Ls;    % Length of squat
g = 9.8;        % Gravity
mu = 0.625;     % Damping Coefficent

datumPE = Ls; % PE Datum
Hs = datumPE-Ls/2;
Hp = datumPE-Lp;

% Moment of inertia
Is = 1/3*Ms*Ls^2;    % Moment of inertia swing
Ip = Mp*Lp.^2;       % Moment of inertia person

omegaN = sqrt(g*(Ms*Ls/2+Mp*Lsq)/(Is+Isq));
mu = 2*Damping*omegaN;

% Initial Conditions
thetaI = thetaN(1);                    % Initial angle (rad)
thetaDotI = thetaDotSN(1);              % Initial velocity (rad/s)

% Time parameters
odeTime = 0.8;                          % Simulation and overflow time for ODE45
nStep = size(storeRun,1);               % Number of total cycles
n = 50;                                 % Number of undamped values per cycle
numCycle = 1;

%% Plot Experimental Diagrams
plotExperimental(dataInput)

%% Initalize Variables
% Compile States
Swing = [Ms Ls Hs Is];

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
if storeRun(find(storeRun(:,3),1),3)==2
    nrgTotal(1) = nrgT(thetaI,thetaDotI,Ip(1),Lp(1),Hp(1));
    stage = 1;
elseif storeRun(find(storeRun(:,3),1),3)==1
    nrgTotal(1) = nrgT(thetaI,thetaDotI,Ip(1),Lp(1),Hp(1));
    stage = 2;
end

index = 0;
tic
for i = 1:2:nStep-2
    index = index+1;
    %% Gaining Kineteic Energy
    if stage == 1 % Gaining kinetic energy need to squat
        
        % Ending System Properties
        endCond = [thetaN(storeRun(i,2)) 7 thetaDotSN(storeRun(i,2))/abs(thetaDotSN(storeRun(i,2)))];
        Squat = [Mp Lp(i) Lp(i) Lp(i)];
        Stand = [Mp Lp(i+2) Lp(i+2) Lp(i+2)];
        moveInput = [mu g moveTime(i+1)];

        % Run ODE
        [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(index,1),thetaDotStep(index,1),moveInput,Swing,[Squat; Squat],odeTime,endCond);   % Stationary Damped

        [thetaMove,thetaDotMove,timeMove] = runODE(thetaDamp(end),thetaDotDamp(end),moveInput,Swing,[Squat; Stand],moveTime(i+1),[0 3]);   % Moving Damped

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
        thetaStep(index+1,:) = real(thetaMove(end));
        thetaDotStep(index+1,:) = real(thetaDotMove(end));
        storeState(end+1:end+size(timeDamp,1)-1,1) = 0; % Stationary
        storeState(end+1:end+size(timeMove,1)-1,1) = 2; % Squatting
        storeTime(end+1:end+size(timeDamp,1)-1,1) = timeDamp(2:end);
        storeTime(end+1:end+size(timeMove,1)-1,1) = timeMove(2:end);
        storeTheta(end+1:end+size(thetaDamp,1)-1,1) = thetaDamp(2:end);
        storeTheta(end+1:end+size(thetaMove,1)-1,1) = thetaMove(2:end);
        storeThetaDot(end+1:end+size(thetaDotDamp,1)-1,1) = thetaDotDamp(2:end);
        storeThetaDot(end+1:end+size(thetaDotMove,1)-1,1) = thetaDotMove(2:end);

        nrgTotal(index+1) = nrgT(thetaMove(end),thetaDotMove(end),Ip(i),Lp(i),Hp(i));
        stage = 2;

    %% Gaining Potential Energy
    elseif stage == 2 % Gaining potential energy need to stand
        
        % Ending System Properties
        endCond = [thetaDotSN(storeRun(i,2)) 8 -thetaN(storeRun(i,2))/abs(thetaN(storeRun(i,2)))];
        Squat = [Mp Lp(i) Hp(i) Ip(i)];
        Stand = [Mp Lp(i+2) Hp(i+2) Ip(i+2)];
        moveInput = [mu g moveTime(i+1)];
        
        % Run ODE
        [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(index,1),thetaDotStep(index,1),moveInput,Swing,[Stand; Stand],odeTime,endCond);   % Stationary Damped

        [thetaMove,thetaDotMove,timeMove] = runODE(thetaDamp(end),thetaDotDamp(end),moveInput,Swing,[Squat; Stand],moveTime(i+1),[0 3]);   % Moving Damped

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
        thetaStep(index+1,:) = real(thetaMove(end));
        thetaDotStep(index+1,:) = real(thetaDotMove(end));
        storeState(end+1:end+size(timeDamp,1)-1,1) = 0; % Stationary
        storeState(end+1:end+size(timeMove,1)-1,1) = 1; % Standing
        storeTime(end+1:end+size(timeDamp,1)-1,1) = timeDamp(2:end);
        storeTime(end+1:end+size(timeMove,1)-1,1) = timeMove(2:end);
        storeTheta(end+1:end+size(thetaDamp,1)-1,1) = thetaDamp(2:end);
        storeTheta(end+1:end+size(thetaMove,1)-1,1) = thetaMove(2:end);
        storeThetaDot(end+1:end+size(thetaDotDamp,1)-1,1) = thetaDotDamp(2:end);
        storeThetaDot(end+1:end+size(thetaDotMove,1)-1,1) = thetaDotMove(2:end);

        nrgTotal(index+1) = nrgT(thetaMove(end),thetaDotMove(end),Ip(i),Lp(i),Hp(i));
        stage = 1;
    end

    nrgCheck = nrgTotal(index+1)-nrgT(0,0,Ip(i),Lp(i),Hp(i));

    if abs(thetaStep(index+1,1))>=pi
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
    elseif index == 150
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
ylim([-1000 1000]*pi/180)

figure(2)
subplot(2,2,1)
title('Angular Displacement')
xlabel('Time (s)')
ylabel('Angular Displacement (rad)')
xlim([0 timeN(end)])
ylim([-4 4])

subplot(2,2,3)
title('Angular Velocity')
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')
xlim([0 timeN(end)])
ylim([-15 15])

%% Create Energy Plots
figure
hold on
plot(nrgTotal(2:2:end),'b')
plot(nrgTotal(1:2:end),'r')
legend('Max Potential Energy','Max Kinetic Energy');
xlabel('Half-Oscillations')

nrgOver = nrgT(pi,0,Ip(end),Lp(end),Hp(end));
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

% figure
% plot(thetaN)
periodSteps = 55;

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

toc





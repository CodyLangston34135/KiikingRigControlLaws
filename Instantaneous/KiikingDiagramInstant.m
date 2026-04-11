%% Kiiking Diagrams
% Cody Langston
clc; clear all; close all;
addpath(fullfile(pwd, 'DiagramProcesses'));

%% User Inputs
% Parameters
Ms = 0.157;     % Mass of swing
Mp = 0.3;       % Mass of person
Ls = 0.305;     % Length of swing
Lsq = 0.305;    % Length of squat
Lst = 0.28;   % Length of stand
g = 9.8;
mu = 1.5;

datumPE = Ls; % PE Datum
Hs = datumPE-Ls/2;
Hsq = datumPE-Lsq;
Hst = datumPE-Lst;

% Moment of inertia
Is = 1/3*Ms*Ls^2;   % Moment of inertia swing
Isq = Mp*Lsq^2;     % Moment of inertia squat
Ist = Mp*Lst^2;     % Moment of inertia stand

% Initial Conditions
thetaI = pi/4;              % Initial angle (rad)
thetaDotI = 0;              % Initial velocity (rad/s)

% Time parameters
odeTime = [0.4,1];          % Simulation and overflow time for ODE45
nStep = 150;                % Number of total cycles
n = 200;                    % Number of undamped values per cycle
moveTime = 1;
numCycle = 1;
% options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3);
options = odeset('MaxStep', 1e-1, 'InitialStep', 1e-1);

% Compile States
Swing = [Ms Ls Hs Is];
Squat = [Mp Lsq Hsq Isq];
Stand = [Mp Lst Hst Ist];

% Initialize Variables
thetaStep = [thetaI thetaI thetaI];                     % [Damped Gain Loss]
thetaDotStep = [thetaDotI thetaDotI thetaDotI];         % [Damped Gain Loss]
nrgStep = [0 0 0];                                      % [PEsq=KEsq=PEi KEst=PEst PEf]
endCond = [0 0];                                        % [endValue endType] 1 = theta, 2 = thetaDot
overflow = 0;                                           % 0 = no overflow, 1 = damps out, 2 = goes over

%% Energy Calculations
nrgKEr = @(thetaDot,I) 1/2*I*thetaDot.^2;               % Rotational Kinetic Energy
nrgPE = @(theta,M,L,H) M*g*(H+L*(1-cos(theta)));        % Gravitational Potential Energy

% Total energy
nrgT = @(theta,thetaDot,Ip,Lp,Hp) (nrgPE(theta,Mp,Lp,Hp)+nrgPE(theta,Ms,Ls/2,Hs))...    % Potential
     + (nrgKEr(thetaDot,Ip)+nrgKEr(thetaDot,Is));                                       % Rotational Kinetic

figure(1)
xlabel('Theta')
ylabel('Theta Dot')
hold on
index = 1;

tic
for i = 1:nStep

    %% Gaining Kineteic Energy
    if thetaDotStep(i,1)==0 || thetaDotStep(i,1)*thetaStep(i,1)<0
       
        % Initial Energy
        nrgPEiSq = nrgT(thetaStep(i,1),thetaDotStep(i,1),Isq,Lsq,Hsq);
        nrgStep(index,1) = nrgPEiSq;

        % Ending System Properties
        endCond = [0 1];
        % [thetaGain,thetaDotGain] = gainKinetic(thetaStep(i,1),thetaDotStep(i,1),endCond(1),g,Swing,Squat,n);                % Undamped Ideal
        % [thetaLoss,thetaDotLoss] = gainKinetic(thetaStep(i,1),thetaDotStep(i,1),endCond(1),g,Swing,Stand,n);                % Undamped Non-Ideal
        [thetaGain,thetaDotGain] = runODE(thetaStep(i,1),thetaDotStep(i,1),endCond,numCycle,0,g,Swing,[Squat; Squat],moveTime,odeTime,options);        % Undamped Ideal
        [thetaLoss,thetaDotLoss] = runODE(thetaStep(i,1),thetaDotStep(i,1),endCond,numCycle,0,g,Swing,[Stand; Stand],moveTime,odeTime,options);        % Undamped Non-Ideal

        [thetaDamp,thetaDotDamp,time,overflow] = runODE(thetaStep(i,1),thetaDotStep(i,1),endCond,numCycle,mu,g,Swing,[Squat; Squat],moveTime,odeTime,options);   % Damped Ideal

        % Phase Diagram
        figure(1)
        plot(thetaGain,thetaDotGain,'--b');
        plot(thetaLoss,thetaDotLoss,'--r')
        plot(thetaDamp,thetaDotDamp,'g');

        % Momentum Gain
        thetaDotGainF = conserveMomentum(thetaDotGain(end),Swing,Squat,Stand);
        thetaDotDampF = conserveMomentum(thetaDotDamp(end),Swing,Squat,Stand);

        % Phase Diagram
        plot([thetaGain(end) thetaGain(end)],[thetaDotGain(end) thetaDotGainF],'--b')
        plot([thetaGain(end) thetaGain(end)],[thetaDotDamp(end) thetaDotDampF],'g')

        % Store Values
        thetaStep(i+1,:) = real([thetaDamp(end) thetaGain(end) thetaLoss(end)]);
        thetaDotStep(i+1,:) = real([thetaDotDampF thetaDotGainF thetaDotLoss(end)]);

    %% Gaining Potential Energy
    elseif thetaStep(i,1)==0 || thetaDotStep(i,1)*thetaStep(i,1)>0 % Gaining kinetic energy

        % Initial Energy
        nrgKEiSt = nrgT(thetaStep(i,2),thetaDotStep(i,2),Ist,Lst,Hst);
        nrgStep(index,2) = nrgKEiSt;

        % Ending System Properties
        endCond = [0 2];
        % [thetaGain,thetaDotGain] = gainPotential(thetaStep(i,2),thetaDotStep(i,2),endCond(1),g,Swing,Stand,n);              % Undamped Ideal
        % [thetaLoss,thetaDotLoss] = gainPotential(thetaStep(i,3),thetaDotStep(i,3),endCond(1),g,Swing,Stand,n);              % Undamped Non-Ideal
        [thetaGain,thetaDotGain] = runODE(thetaStep(i,2),thetaDotStep(i,2),endCond,numCycle,0,g,Swing,[Stand; Stand],moveTime,odeTime,options);        % Undamped Ideal
        [thetaLoss,thetaDotLoss] = runODE(thetaStep(i,3),thetaDotStep(i,3),endCond,numCycle,0,g,Swing,[Stand; Stand],moveTime,odeTime,options);        % Undamped Non-Ideal

        [thetaDamp,thetaDotDamp,time,overflow] = runODE(thetaStep(i,1),thetaDotStep(i,1),endCond,numCycle,mu,g,Swing,[Stand; Stand],moveTime,odeTime,options);   % Damped Ideal

        % Phase Diagram
        plot(thetaGain,thetaDotGain,'--b');
        plot(thetaLoss,thetaDotLoss,'--r')
        plot(thetaDamp,thetaDotDamp,'g');

        % Final Damped Energy
        nrgPEd = nrgT(thetaDamp(end),thetaDotDamp(end),Ist,Lst,Hst);
        nrgStep(index,3) = nrgPEd;
        index = index + 1;
        
        % Store Values
        thetaStep(i+1,:) = real([thetaDamp(end) thetaGain(end) thetaLoss(end)]);
        thetaDotStep(i+1,:) = real([thetaDotDamp(end) thetaDotGain(end) thetaDotLoss(end)]);

    end

    if thetaStep(i+1,1) == 0 && thetaDotStep(i+1,1) == 0 || overflow == 1
        disp('Damps Out')
        overflow = 1;
        nrgStep(end,:) = [];
        break
    elseif abs(thetaStep(i+1,1))>=pi || overflow == 2
        disp('Goes Over')
        overflow = 2;

        xline(pi,'--k')
        xline(-pi,'--k')
        xlim([-pi-0.5 pi+0.5])

        % Final Kinetic Energy
        nrgKEf = nrgT(0,max(abs(thetaDotDamp)),Ist,Lst,Hst);
        nrgStep(end,3) = nrgKEf;
        break
    end
end

nrgOver = nrgT(pi,0,Ist,Lst,Hst);
% nrgStep(end,:) = [];
nrgGain = nrgStep(:,2)-nrgStep(:,1);
nrgLost = nrgStep(:,2)-nrgStep(:,3);
nrgTotal = nrgStep(:,3);

figure(2)
hold on
plot(nrgGain,'--b')
plot(nrgLost,'--r')
plot(nrgTotal,'g')
if overflow == 2
    yline(nrgOver,'--k');
    legend('Energy Gained','Energy Lost','System Energy','Required Energy')
else
    legend('Energy Gained','Energy Lost','System Energy')
end


toc
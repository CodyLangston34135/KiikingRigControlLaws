%% Kiiking Model
% Cody Langston
clc; close all; clear all;

%% User Input
% Parameters
Ms = 0.157;     % Mass of swing
Mp = 0.328;     % Mass of person
Ls = 0.305/2;   % Length of swing
Lsq = 0.28;     % Length of squat
Lst = 0.25;     % Length of stand
g = 9.8;
% mu = 0.000225;
mu = 0.0;

datumPE = Ls*2; % PE Datum
Hs = datumPE-Ls;
Hsq = datumPE-Lsq;
Hst = datumPE-Lst;

% Moment of inertia
Is = 1/3*Ms*Ls^2;   % Moment of inertia swing
Isq = Mp*Lsq^2;     % Moment of inertia squat
Ist = Mp*Lst^2;     % Moment of inertia stand

% Initial Conditions
thetaI = pi/4;            % Initial angle (rad)
thetaDotI = 0;            % Initial velocity (rad/s)
statePos = [0 0 0 0 0];     % [theta thetaDot nrgGain state]

% Time parameters
odeTime = 0.5;                % Simulation time for ODE45
nStep = 6;                 % Number of quarter cycles
n = 100;                    % Number of undamped values per cycle

%% Energy Calculations
nrgKEr = @(thetaDot,I) 1/2*I*thetaDot.^2;               % Rotational Kinetic Energy
nrgKEt = @(thetaDot,M,L) 1/2*M*(thetaDot*L).^2;         % Translational Kinetic Energy
nrgPE = @(theta,M,L,H) M*g*(H+L*(1-cos(theta)));        % Gravitational Potential Energy

% Total energy
nrgT = @(theta,thetaDot,Ip,Lp,Hp) (nrgPE(theta,Mp,Lp,Hp)+nrgPE(theta,Ms,Ls,Hs))... % Potential
     + (nrgKEr(thetaDot,Ip)+nrgKEr(thetaDot,Is))...                                % Rotational Kinetic
     + (nrgKEt(thetaDot,Mp,Lp)+nrgKEt(thetaDot,Ms,Ls));                            % Translational Kinetic

calcTheta = @(nrgT,nrgKEf,Lp,Hp) acos(1-((nrgT-nrgKEf)-(Ms*g*Hs+Mp*g*Hp))/(Ms*g*Ls+Mp*g*Lp));
calcThetaDot = @(nrgT,nrgPEf,Lp,Ip) sqrt((nrgT-nrgPEf)/(1/2*(Is+Ip+Ms*Ls^2+Mp*Lp^2)));

%% Initialize variables
statePos = [0 0 0 0 0];

%% Calculations
% Gaining kinetic or potential energy
if thetaI==0 || thetaDotI*thetaI>0 % Gaining potential energy
    state = 1;
    nrgTI = nrgT(thetaI,thetaDotI,Isq,Lsq,Hsq);
    storeCycle = [thetaI thetaDotI thetaI thetaDotI thetaI thetaDotI];
elseif thetaDotI==0 || thetaDotI*thetaI<0 % Gaining kinetic energy
    state = 2;
    nrgTI = nrgT(thetaI,thetaDotI,Ist,Lst,Hst);
    storeCycle = [thetaI thetaDotI thetaI thetaDotI thetaI thetaDotI];
end

% Create plot
figure(1)
xlabel('Theta')
ylabel('Theta Dot')
xline(pi,'--k')
hold on
xline(-pi,'--k')

% Store initial energy
storeNRG(1) = nrgTI;
index = 1;

tic
for i = 1:nStep
    % If gaining potential energy
    if state == 1
        % Store previous conditions
        thetaDotL = storeCycle(i,6);
        thetaDotSt = storeCycle(i,4);
        thetaDotDamp = storeCycle(i,2);
        
        thetaDotLN = linspace(0,thetaDotL,n);
        thetaDotStN = linspace(0,thetaDotSt,n);
        thetaDotDampN = linspace(0,thetaDotDamp,n);

        % State 1 Total Energy
        nrgTL = nrgT(0,thetaDotL,Ist,Lst,Hst);
        nrgTst = nrgT(0,thetaDotSt,Ist,Lst,Hst);

        % State 2 Kinetic Energy
        nrgKEL = ((nrgKEr(thetaDotLN,Ist)+nrgKEr(thetaDotLN,Is))...        % Rotational Kinetic
               + (nrgKEt(thetaDotLN,Mp,Lst)+nrgKEt(thetaDotLN,Ms,Ls)));    % Translational Kinetic
        nrgKEst = ((nrgKEr(thetaDotStN,Ist)+nrgKEr(thetaDotStN,Is))...     % Rotational Kinetic
                + (nrgKEt(thetaDotStN,Mp,Lst)+nrgKEt(thetaDotStN,Ms,Ls))); % Translational Kinetic

        % State 2 Anglular Displacement
        thetaL = calcTheta(nrgTL,nrgKEL,Lst,Hst);
        thetaSt = calcTheta(nrgTst,nrgKEst,Lst,Hst);

        % Account for direction of motion
        if thetaDotSt<=0
            thetaL = -thetaL;
            thetaSt = -thetaSt;
        end

        % Damped Method
        [thetaDamp, thetaDotDamp, time] = DampedKiiking([0, thetaDotSt],mu,Ms,Ls,Mp,Lst,g,odeTime);
        thetaDampN = interp1(thetaDotDamp,thetaDamp,thetaDotDampN);
        timeN = interp1(thetaDotDamp,time,thetaDotDampN);
        timeEnd = timeN(end);

        % Plot phase diagram from state 1 to 2
        plot(real(thetaSt),real(thetaDotStN),'--b')
        plot(real(thetaL),real(thetaDotLN),'--r')
        plot(real(thetaDampN),real(thetaDotDampN),'m')

        % Store ending state
        state = 2;
        nrgTstF = nrgT(thetaSt(1),0,Ist,Lst,Hst);
        nrgTDampF = nrgT(thetaDampN(1),0,Ist,Lst,Hst);
        nrgGain = nrgTstF-storeNRG(index);
        nrgLost = nrgTstF-nrgTDampF;
        statePos(i+1,:) = [nrgGain nrgLost nrgTDampF timeEnd state];
        storeCycle(i+1,:) = [thetaDampN(1) 0 thetaSt(1) 0 thetaL(1) 0];

    % If gaining kinetic energy
    elseif state == 2
        % Store previous conditions
        theta = storeCycle(i,1);

        thetaSqN = linspace(0,theta,n);
        thetaStN = linspace(0,theta,n);
        thetaDampN = linspace(0,theta,n);

        % State 2 Total Energy
        nrgSq = nrgPE(theta,Mp,Lsq,Hsq)+nrgPE(theta,Ms,Ls,Hs);
        nrgSt = nrgPE(theta,Mp,Lst,Hst)+nrgPE(theta,Ms,Ls,Hs);

        % State 1 Potential Energy
        nrgPEsq = nrgPE(thetaSqN,Mp,Lsq,Hsq)+nrgPE(thetaSqN,Ms,Ls,Hs);
        nrgPEst = nrgPE(thetaStN,Mp,Lst,Hst)+nrgPE(thetaStN,Ms,Ls,Hs);

        % State 1 Angular Velocity
        thetaDotSq = calcThetaDot(nrgSq,nrgPEsq,Lsq,Isq);
        thetaDotSt = calcThetaDot(nrgSt,nrgPEst,Lst,Ist);

        % Velocity Direction
        if theta>=0
            thetaDotSt = -thetaDotSt;
            thetaDotSq = -thetaDotSq;
        end

        % Damped Method
        [thetaDamp, thetaDotDamp, time] = DampedKiiking([theta, 0],mu,Ms,Ls,Mp,Lsq,g,odeTime);
        thetaDotDampN = interp1(thetaDamp,thetaDotDamp,thetaDampN);
        timeN = interp1(thetaDamp,time,thetaDampN);
        timeEnd = timeN(end);

        % Phase Diagrams
        plot(real(thetaSqN),real(thetaDotSq),'--b')
        plot(real(thetaStN),real(thetaDotSt),'--r')
        plot(real(thetaDampN),real(thetaDotDampN),'m')

        % Store ending state
        state = 3;
        nrgTsqF = nrgT(0,thetaDotSq(1),Isq,Lsq,Hsq);
        nrgTstF = nrgT(0,thetaDotSt(1),Ist,Lst,Hst);
        nrgTDampF = nrgT(0,thetaDotDamp(1),Isq,Lsq,Hsq);
        nrgGain = nrgTsqF-nrgTstF;
        nrgLost = nrgTsqF-nrgTDampF;
        statePos(i+1,:) = [nrgGain nrgLost nrgTDampF timeEnd state];

        storeCycle(i+1,:) = [0 thetaDotDampN(1) 0 thetaDotSq(1) 0 thetaDotSt(1)];

        % Store energy for entire cycle
        index = index+1;
        storeNRG(index) = nrgTstF;
   
    % If changing position
    elseif state ==3
        % Store previous condition
        thetaDotSq = storeCycle(i,4);
        thetaDotDampSq = storeCycle(i,2);
        timeEnd = statePos(i,4);
        
        % Conservation of Angular Momentum
        thetaDotSt = thetaDotSq*(Isq+Is)/(Ist+Is);
        thetaDotDampSt = thetaDotDampSq*(Isq+Is)/(Ist+Is);

        % Plot change due to angular momentum
        plot([0 0], real([thetaDotSq, thetaDotSt]),'--b')
        plot([0 0], real([thetaDotDampSq, thetaDotDampSt]),'m')

        % Store ending state
        state = 1;
        nrgTsqF = nrgT(0,thetaDotSq,Isq,Lsq,Hsq);
        nrgTstF = nrgT(0,thetaDotSt,Ist,Lst,Hst);
        nrgTDampF = nrgT(0,thetaDotDampSt,Ist,Lst,Hst);
        nrgGain = nrgTstF-storeNRG(index);
        nrgLost = nrgTstF-nrgTDampF;

        statePos(i+1,:) = [nrgGain nrgLost nrgTDampF timeEnd state];
        storeCycle(i+1,:) = [0 thetaDotDampSt 0 thetaDotSt 0 storeCycle(i,6)];
    else
        disp('Non-zero starting position')
    end
end

toc

statePos;


function [theta, thetaDot, time] = DampedKiiking(IC,mu,Ms,Ls,Mp,Lp,g,tFinal)

    I = 1/3*Ms*Ls^2+Mp*Lp^2;
    ML = Mp*Lp + Ms*Ls/2;
    
    thDot = @(t,th) [th(2); (-mu*sign(th(2))-ML*g*sin(th(1))/I)];
    % thDot = @(t,th) [th(2); (-mu*(th(2))-ML*g*sin(th(1))/I)];
    tSpan = [0 tFinal];
    
    options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3);
    
    [time, thetaOut] = ode45(thDot,tSpan, IC);
    
    theta = thetaOut(:,1);
    thetaDot = thetaOut(:,2);

end
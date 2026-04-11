%% Kiiking Model
% Cody Langston
clc; close all; clear all;

%% User Input
% Parameters
Ms = 0.157;     % Mass of swing
Mp = 0.328;     % Mass of person
Ls = 0.305/2;   % Length of swing
Lsq = 0.305;     % Length of squat
Lst = 0.28;     % Length of stand
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
nStep = 20;                 % Number of quarter cycles
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

%% Calculations
% Gaining kinetic or potential energy
if thetaI==0 || thetaDotI*thetaI>0 % Gaining potential energy
    state = 1;
    nrgTI = nrgT(thetaI,thetaDotI,Isq,Lsq,Hsq);
    statePos(1,:) = [thetaI thetaDotI 0 0 state];
elseif thetaDotI==0 || thetaDotI*thetaI<0 % Gaining kinetic energy
    state = 2;
    nrgTI = nrgT(thetaI,thetaDotI,Ist,Lst,Hst);
    statePos(1,:) = [thetaI thetaDotI 0 0 state];
end

% Create plot
figure(1)
xlabel('Theta')
ylabel('Theta Dot')
xline(pi,'--k')
hold on
xline(-pi,'--k')

tic
for i = 1:nStep
    % If gaining potential energy
    if state == 1
        % Store previous conditions
        thetaDotSq = statePos(i-1,2);
        thetaDotSt = statePos(i,2);
        
        thetaDotSqN = linspace(0,thetaDotSq,n);
        thetaDotStN = linspace(0,thetaDotSt,n);
        thetaDotDampN = linspace(0,thetaDotSt,n);

        % State 1 Total Energy
        nrgTsq = nrgT(0,thetaDotSq,Isq,Lsq,Hsq);
        nrgTst = nrgT(0,thetaDotSt,Ist,Lst,Hst);

        nrgKEst2 = ((nrgKEr(thetaDotSt,Ist)+nrgKEr(thetaDotSt,Is))...     % Rotational Kinetic
                + (nrgKEt(thetaDotSt,Mp,Lst)+nrgKEt(thetaDotSt,Ms,Ls))); % Translational Kinetic

        nrgPEst2 = (nrgPE(0,Mp,Lst,Hst)+nrgPE(0,Ms,Ls,Hs))


        % State 2 Kinetic Energy
        nrgKEsq = ((nrgKEr(thetaDotSqN,Isq)+nrgKEr(thetaDotSqN,Is))...     % Rotational Kinetic
                + (nrgKEt(thetaDotSqN,Mp,Lsq)+nrgKEt(thetaDotSqN,Ms,Ls))); % Translational Kinetic
        nrgKEst = ((nrgKEr(thetaDotStN,Ist)+nrgKEr(thetaDotStN,Is))...     % Rotational Kinetic
                + (nrgKEt(thetaDotStN,Mp,Lst)+nrgKEt(thetaDotStN,Ms,Ls))); % Translational Kinetic

        % State 2 Anglular Displacement
        thetaSq = calcTheta(nrgTsq,nrgKEsq,Lsq,Hsq);
        thetaSt = calcTheta(nrgTst,nrgKEst,Lst,Hst);

        % Account for direction of motion
        if thetaDotSt<=0
            thetaSq = -thetaSq;
            thetaSt = -thetaSt;
        end

        % Damped Method
        [thetaDamp, thetaDotDamp] = DampedKiiking([0, thetaDotSt],mu,Ms,Ls,Mp,Lst,odeTime);
        thetaDampN = interp1(thetaDotDamp,thetaDamp,thetaDotDampN);

        % Plot phase diagram from state 1 to 2
        plot(real(thetaSt),real(thetaDotStN),'--r')
        plot(real(thetaSq),real(thetaDotSqN),'--b')
        plot(real(thetaDampN),real(thetaDotDampN),'m')

        % Store ending state
        state = 2;
        nrgTstF = nrgT(thetaSt(1),0,Ist,Lst,Hst);
        nrgTsqF = nrgT(thetaSq(1),0,Isq,Lsq,Hsq);
        nrgTDampF = nrgT(thetaDampN(1),0,Ist,Lst,Hst);
        nrgGain = nrgTstF-nrgTsqF;
        nrgLost = nrgTstF-nrgTDampF;
        statePos(i+1,:) = [thetaDampN(1) 0 nrgGain nrgLost state];

    % If gaining kinetic energy
    elseif state == 2
        % Store previous conditions
        theta = statePos(i,1);

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
        [thetaDamp, thetaDotDamp] = DampedKiiking([theta, 0],mu,Ms,Ls,Mp,Lsq,odeTime);
        thetaDotDampN = interp1(thetaDamp,thetaDotDamp,thetaDampN);

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
        statePos(i+1,:) = [0 thetaDotDampN(1) nrgGain nrgLost state];
        storeState3 = thetaDotSq(1);
   
    % If changing position
    elseif state ==3
        % Store previous condition
        thetaDotSq = storeState3;
        thetaDotDampSq = statePos(i,2);
        
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
        nrgGain = nrgTstF-nrgTsqF;
        nrgLost = nrgTstF-nrgTsqF;

        statePos(i+1,:) = [0 thetaDotDampSt nrgGain nrgLost state];
    else
        disp('Non-zero starting position')
    end
end

toc
statePos;


function [theta, thetaDot, time] = DampedKiiking(IC,mu,Ms,Ls,Mp,Lp,tFinal)

    g = 9.8;
    A = 1/3*Ms*Ls^2+Mp*Lp^2;
    B = Mp*Lp + Mp*Lp/2;
    
    % thDot = @(t,th) [th(2); (-mu*sign(th(2))-B*sin(th(1))/A)];
    thDot = @(t,th) [th(2); (-mu*(th(2))-B*g*sin(th(1))/A)];
    tSpan = [0 tFinal];
    
    % options = odeset('MaxStep', 1e-3, 'InitialStep', 1e-3);
    
    [time, thetaOut] = ode45(thDot,tSpan, IC);
    
    theta = thetaOut(:,1);
    thetaDot = thetaOut(:,2);

end
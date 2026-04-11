%% Kiiking Analytical Diagrams
% Cody Langston
clc; clear all; %close all;
addpath(fullfile(pwd, 'NumericalProcesses'));

%% User Inputs

Mass = linspace(0.45,0,400);    % High mass to low mass
% Length = linspace(0.2,100,200);
Damping = linspace(0.1,0,400); % High damping to low damping
% Damping = linspace(3,0,200);

debug = false;  % Debug energy levels with ODEs

for i = 1:length(Damping)   % Loop damping        
    for j = 1:length(Mass)  % Loop mass
        % Parameters
        Ms = 0.157;         % Mass of swing
        Mp = Mass(j);       % Mass of person
        % Mp = 0.3;
        % Ls = Length(j);
        Ls = 0.305;         % Length of swing
        Lsq = Ls;           % Length of squat
        Lst = Ls-0.025;     % Length of stand
        g = 9.8;
        
        % Moment of inertia
        Is = 1/3*Ms*Ls^2;   % Moment of inertia swing
        Isq = Mp*Lsq^2;     % Moment of inertia squat
        Ist = Mp*Lst^2;     % Moment of inertia stand
        
        % Combinations
        MLsq = (Ms*Ls/2+Mp*Lsq);
        MLst = (Ms*Ls/2+Mp*Lst);
        MIsq = Is+Isq;
        MIst = Is+Ist;

        % Potential Energy Datum
        datumPE = Ls;
        Hs = datumPE-Ls/2;
        Hsq = datumPE-Lsq;
        Hst = datumPE-Lst;

        % Modal Parameters
        omegaN = sqrt(g*MLsq/MIsq); % Undamped natural frequency
        mu = 2*Damping(i)*omegaN;   % Damping coefficient
        
        % Compile States
        Swing = [Ms Ls Hs Is];
        Squat = [Mp Lsq Hsq Isq];
        Stand = [Mp Lst Hst Ist];
        
        % Initialize Variables
        theta = zeros(0,1);
        thetaDot = zeros(0,1);
        
        %% Energy Calculations
        nrgKE = @(thetaDot,I) 1/2*I*thetaDot.^2;               % Rotational Kinetic Energy
        nrgPE = @(theta,M,L) M*g*(L*(1-cos(theta)));        % Gravitational Potential Energy

        % Viscous Damping
        nrgDampPE = @(theta,x,Mp,mu) -mu*trapz(x,sqrt(2*g*MLsq/MIsq*(cos(x)-cos(theta))));
        nrgDampKE = @(thetaDot,x,Mp,mu) -mu*trapz(x,sqrt((1/2*MIst*thetaDot^2-g*MLst*(1-cos(x)))/(1/2*MIst)));
        nrgDampPEcumul = @(theta,x,Mp,mu) -mu*cumtrapz(x,sqrt(2*g*MLsq/MIsq*(cos(x)-cos(theta))));
        nrgDampKEcumul = @(thetaDot,x,Mp,mu) -mu*cumtrapz(x,sqrt((1/2*MIst*thetaDot^2-g*MLst*(1-cos(x)))/(1/2*MIst)));

        % Coulomb Damping
        % nrgDampPE = @(theta,x,Mp,mu) -mu*trapz(x,sign(sqrt(2*g*MLsq/MIsq*(cos(x)-cos(theta)))));
        % nrgDampKE = @(thetaDot,x,Mp,mu) -mu*trapz(x,sign(sqrt((1/2*MIst*thetaDot^2-g*MLst*(1-cos(x)))/(1/2*MIst))));
        % nrgDampPEcumul = @(theta,x,Mp,mu) -mu*cumtrapz(x,sign(sqrt(2*g*MLsq/MIsq*(cos(x)-cos(theta)))));
        % nrgDampKEcumul = @(thetaDot,x,Mp,mu) -mu*cumtrapz(x,sign(sqrt((1/2*MIst*thetaDot^2-g*MLst*(1-cos(x)))/(1/2*MIst))));

        % Total energy
        nrgT = @(theta,thetaDot,Ip,Lp) (nrgPE(theta,Mp,Lp)+nrgPE(theta,Ms,Ls/2))...         % Potential
             + (nrgKE(thetaDot,Ip)+nrgKE(thetaDot,Is));                                     % Rotational Kinetic
        
        % Create Analytical Expressions
        calcThetaDot2 = @(theta1,Mp,mu,x) sqrt((nrgPE(theta1,Ms,Ls/2)+nrgPE(theta1,Mp,Lsq)-abs(MIsq*nrgDampPE(theta1,x,Mp,mu)))/(1/2*MIsq));
        calcThetaDot3 = @(thetaDot2,Mp,mu) thetaDot2*MIsq/MIst;
        calcNRG4 = @(thetaDot3,thetaEnd,Mp,mu,x) (nrgKE(thetaDot3,Is)+nrgKE(thetaDot3,Ist)-abs(MIst*nrgDampKE(thetaDot3,x,Mp,mu)));
        calcNRGEnd = @(theta1,Mp,mu) calcNRG4(calcThetaDot3(calcThetaDot2(theta1)),theta1);

        %% Calculate if swing gains energy at pi
        theta1 = pi;
        x = linspace(0,theta1,100);
        nrgOver = nrgPE(theta1,Ms,Ls/2)+nrgPE(theta1,Mp,Lst);
        thetaDot2 = calcThetaDot2(theta1,Mp,mu,x);
        thetaDot3 = calcThetaDot3(thetaDot2,Mp,mu);
        nrg4 = calcNRG4(thetaDot3,pi,Mp,mu,x);

        % Debug States
        if debug
            thetaN = linspace(theta1,0,100);
            thetaDotUndamped = sqrt(2*g*(Ms*Ls/2+Mp*Lsq)/(1/3*Ms*Ls^2+Mp*Lsq^2)*(cos(thetaN)-cos(theta1)));
            workDamp = MIsq*nrgDampPEcumul(theta1,thetaN,Mp,mu);
            thetaDotDamped = sqrt(2*(nrgPE(theta1,Ms,Ls/2)+nrgPE(theta1,Mp,Lsq)-nrgPE(thetaN,Ms,Ls/2)-nrgPE(thetaN,Mp,Lsq)-workDamp)/(MIsq));
    
            endCond = [0 1];
            moveTime = 0.1;
            odeTime = 2;
            moveInput = [mu g moveTime];
            [thetaDamp,thetaDotDamp,timeDamp] = runODE(pi-0.01,0,moveInput,Swing,[Squat; Squat],odeTime,endCond);   % Damped Ideal
            thetaDotDampN = interp1(thetaDamp,thetaDotDamp,thetaN);
            thetaDotDampN(1) = 0;
    
            nrgDamp = nrgKE(thetaDotDampN,Isq)-nrgKE(thetaDotUndamped,Isq);
    
            plot(thetaN,-thetaDotUndamped)
            hold on
            plot(thetaN,-thetaDotDamped)
            plot(thetaN,thetaDotDampN)

            legend('Undamped','Analytical Damped','ODE45 Damped')

            figure
            plot(thetaN,workDamp)
            hold on
            plot(thetaN,-nrgDamp)
            plot(thetaN,workDamp./nrgDamp)

            thetaDotN3 = thetaDotDamped(end)*MIsq/MIst;

            thetaN = linspace(0,theta1,100);
            thetaDotUndamped = sqrt((1/2*MIst*thetaDotN3^2-g*MLst*(1-cos(thetaN)))/(1/2*MIst));
           
            % workDamp = -sign(mu)*MIst*cumtrapz(thetaN,thetaDotUndamped);
            workDamp = MIst*nrgDampKEcumul(thetaDotN3,thetaN,Mp,mu);
            thetaDotDamped = sqrt((nrgKE(thetaDotN3,MIst)-nrgPE(thetaN,Ms,Ls/2)-nrgPE(thetaN,Mp,Lst)+workDamp)/(1/2*MIst));
    
            endCond = [0 2];
            moveTime = 0.1;
            odeTime = 2;
            moveInput = [mu g moveTime];
            [thetaDamp,thetaDotDamp,timeDamp] = runODE(0,thetaDotN3,moveInput,Swing,[Stand; Stand],odeTime,endCond);   % Damped Ideal
            thetaDotDampN = interp1(thetaDamp,thetaDotDamp,thetaN);
    
            nrgDamp = nrgKE(thetaDotDampN,Ist)-nrgKE(thetaDotUndamped,Ist);
            nrgTotal = nrgT(thetaN,thetaDotDampN,Ist,Lst);
    
            figure
            plot(-thetaN,-thetaDotUndamped)
            hold on
            plot(-thetaN,-thetaDotDamped)
            plot(-thetaN,-thetaDotDampN)

            legend('Undamped','Analytical Damped','ODE45 Damped')

            figure
            plot(thetaN,-workDamp)
            hold on
            plot(thetaN,-nrgDamp)
            plot(thetaN,workDamp./nrgDamp)
            
        end

        if nrg4>=nrgOver && thetaDot2 >= 0 && thetaDot3 >= 0 % If swing gains energy store going over
            % disp('Goes over')
            state = pi;
        else % If swing loses energy dont store going over
            error = 1;
            theta = [pi pi/16];  % Minimize change in energy between pi/4 and pi
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
        
            if thetaTemp <= theta(2)+.2     % if LCO is near lower bound, damps out
                % disp('Damps Out')
                state = thetaTemp;
            else                            % if LCO was found, store amplitude
                % disp('Limit Cycle Oscillation')
                state = thetaTemp;
            end
        end
        storeState(j,i) = state;
    end
end
figure
surf(Damping,Mass,storeState)
shading interp
% colormap jet

xlabel('Damping');
ylabel('Mass');
zlabel('Theta')
title('Analytical')




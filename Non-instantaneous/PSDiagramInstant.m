%% Kiiking Optimized Diagrams
% Cody Langston
clc; clear all; close all;
addpath(fullfile(pwd, 'NumericalProcesses'));
addpath(fullfile(pwd, 'Figures'));
% pool = parpool("Threads");

Mass = linspace(2,0,40);
% Length = linspace(0.1,5,10);

Damping = linspace(0,0.1,40); % Viscous Damping
dampTheta = pi/6;             % Viscous Damping
% Damping = linspace(0,3,40);    % Coulomb Damping
% dampTheta = pi/4;              %Coulomb Damping
plotFigure = false;
killDampedOut = true;

storeOverflow = ones(length(Mass),length(Damping));

tic
for k = 1:length(Damping)
    fprintf('Damping %.2f\n',Damping(k))
    for l = 1:length(Mass)
        if storeOverflow(l,k) == 0
        else
        %% User Inputs
        % Parameters
        Ms = 0.157;     % Mass of swing
        Mp = Mass(l);       % Mass of person
        % Mp = 0.3;
        % Ls = Length(j);
        Ls = 0.305;     % Length of swing
        Lsq = Ls;    % Length of squat
        Lst = Ls-0.025;   % Length of stand
        g = 9.8;

        % Mp = 0.3;
        % Ls = Length(l);
        % Lsq = Ls;
        % Lst = Ls-0.015;
        
        datumPE = Ls; % PE Datum
        Hs = datumPE-Ls/2;
        Hsq = datumPE-Lsq;
        Hst = datumPE-Lst;
        
        % Moment of inertia
        Is = 1/3*Ms*Ls^2;   % Moment of inertia Swing
        Isq = Mp*Lsq^2;     % Moment of inertia Squat
        Ist = Mp*Lst^2;     % Moment of inertia Stand

        omegaN = sqrt(g*(Ms*Ls/2+Mp*Lsq)/(Is+Isq));
        mu = 2*Damping(k)*omegaN;
        
        % Initial Conditions
        thetaI = pi/2;              % Initial angle (rad)
        thetaDotI = 0;              % Initial velocity (rad/s)
        
        % Time parameters
        odeTime = 0.8;          % Simulation and overflow time for ODE45
        nStep = 70;             % Number of total cycles
        n = 20;                    % Number of undamped values per cycle
        moveTime = 0.001;
        numCycle = 1;
        
        % Compile States
        Swing = [Ms Ls Hs Is];
        Squat = [Mp Lsq Hsq Isq];
        Stand = [Mp Lst Hst Ist];
        
        % Initialize Variables
        nrgStep = [0 0 0];                                      % [PEsq=KEsq=PEi KEst=PEst PEf]
        endCond = [0 0];                                        % [endValue endType] 1 = theta, 2 = thetaDot
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
        
        if thetaDotI==0 || thetaDotI*thetaI<0
            nrgTotal = nrgT(thetaI,thetaDotI,Isq,Lsq,Hsq);
            stage = 1;
        elseif thetaI ==0 || thetaDotI*thetaI>0
            nrgTotal = nrgT(thetaI,thetaDotI,Ist,Lst,Hst);
            stage = 2;
        end
        
        %% Going Over
        thetaStep = 5*pi/6;
        thetaDotStep = 0;
        for i = 1:7
        
            %% Gaining Kineteic Energy
            if stage == 1 % Gaining kinetic energy need to Squat
        
                % Ending System Properties
                endCond = [0 1];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Squat; Squat],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                % Adjust Values
                timeDamp = timeN(1:end)+currentTime;
                currentTime = timeDamp(end);
                timeN(end+1) = timeN(end);
                thetaDampN(end+1) = thetaDampN(end);
                [thetaDotInstant] = conserveMomentum(thetaDotDampN(end),Swing,Squat,Stand);
                thetaDotDampN(end+1) = thetaDotInstant;

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:end),thetaDotDampN(1:end),'g')
                end

                % Store Values
                thetaStep(i+1,:) = real(thetaDampN(end));
                thetaDotStep(i+1,:) = real(thetaDotDampN(end));
                storeState(end+1:end+size(timeDamp',1)-1,1) = 0; % Stationary
                storeTime(end+1:end+size(timeDamp',1)-1,1) = timeDamp(2:end)';
                storeTheta(end+1:end+length(thetaDampN)-1,1) = thetaDampN(2:end);
                storeThetaDot(end+1:end+length(thetaDotDampN)-1,1) = thetaDotDampN(2:end);

                nrgTotal(i+1) = nrgT(thetaDampN(end),thetaDotDampN(end),Ist,Lst,Hst);
                nrgCheck = nrgTotal(i+1) - nrgT(0,0,Ist,Lst,Hst);
                stage = 2;
        
            %% Gaining Potential Energy
            elseif stage == 2 % Gaining potential energy need to Stand
                
                % Ending System Properties
                endCond = [0 2];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Stand; Stand],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                % Adjust time
                timeDamp = timeN(1:end)+currentTime;
                currentTime = timeDamp(end);

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:end),thetaDotDampN(1:end),'g')
                end
        
                % Store Values
                thetaStep(i+1,:) = real(thetaDampN(end));
                thetaDotStep(i+1,:) = real(thetaDotDampN(end));
                storeState(end+1:end+size(timeDamp',1)-1,1) = 0; % Stationary
                storeTime(end+1:end+size(timeDamp',1)-1,1) = timeDamp(2:end)';
                storeTheta(end+1:end+length(thetaDampN)-1,1) = thetaDampN(2:end);
                storeThetaDot(end+1:end+length(thetaDotDampN)-1,1) = thetaDotDampN(2:end);
        
                nrgTotal(i+1) = nrgT(thetaDampN(end),thetaDotDampN(end),Isq,Lsq,Hsq);
                nrgCheck = nrgTotal(i+1)-nrgT(0,0,Isq,Lsq,Hsq);
                stage = 1;
            end
        end
            
        if nrgTotal(7) >= nrgTotal(5)
            disp('Goes Over')
            overflow = pi;
    
            if plotFigure
                figure(1)
                xline(pi,'--k')
                xline(-pi,'--k')
                xlim([-pi-0.5 pi+0.5])
            end
        else

        %% Not Going Over
        thetaStep = dampTheta;
        thetaDotStep = 0;
        for i = 1:7
        
            %% Gaining Kineteic Energy
            if stage == 1 % Gaining kinetic energy need to Squat
        
                % Ending System Properties
                endCond = [0 1];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Squat; Squat],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                % Adjust Values
                timeDamp = timeN(1:end)+currentTime;
                currentTime = timeDamp(end);
                timeN(end+1) = timeN(end);
                thetaDampN(end+1) = thetaDampN(end);
                [thetaDotInstant] = conserveMomentum(thetaDotDampN(end),Swing,Squat,Stand);
                thetaDotDampN(end+1) = thetaDotInstant;

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:end),thetaDotDampN(1:end),'g')
                end

                % Store Values
                thetaStep(i+1,:) = real(thetaDampN(end));
                thetaDotStep(i+1,:) = real(thetaDotDampN(end));
                storeState(end+1:end+size(timeDamp',1)-1,1) = 0; % Stationary
                storeTime(end+1:end+size(timeDamp',1)-1,1) = timeDamp(2:end)';
                storeTheta(end+1:end+length(thetaDampN)-1,1) = thetaDampN(2:end);
                storeThetaDot(end+1:end+length(thetaDotDampN)-1,1) = thetaDotDampN(2:end);
        
                nrgTotal(i+1) = nrgT(thetaDampN(end),thetaDotDampN(end),Ist,Lst,Hst);
                nrgCheck = nrgTotal(i+1) - nrgT(0,0,Ist,Lst,Hst);
                stage = 2;
        
            %% Gaining Potential Energy
            elseif stage == 2 % Gaining potential energy need to Stand
                
                % Ending System Properties
                endCond = [0 2];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Stand; Stand],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                % Adjust time
                timeDamp = timeN(1:end)+currentTime;
                currentTime = timeDamp(end);

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:end),thetaDotDampN(1:end),'g')
                end
        
                % Store Values
                thetaStep(i+1,:) = real(thetaDampN(end));
                thetaDotStep(i+1,:) = real(thetaDotDampN(end));
                storeState(end+1:end+size(timeDamp',1)-1,1) = 0; % Stationary
                storeTime(end+1:end+size(timeDamp',1)-1,1) = timeDamp(2:end)';
                storeTheta(end+1:end+length(thetaDampN)-1,1) = thetaDampN(2:end);
                storeThetaDot(end+1:end+length(thetaDotDampN)-1,1) = thetaDotDampN(2:end);
        
                nrgTotal(i+1) = nrgT(thetaDampN(end),thetaDotDampN(end),Isq,Lsq,Hsq);
                nrgCheck = nrgTotal(i+1)-nrgT(0,0,Isq,Lsq,Hsq);
                stage = 1;
            end
        end

        if nrgTotal(8) <= nrgTotal(4)
            disp('Damps Out')
            overflow = 0;

            if killDampedOut
                for m = l:length(Mass)
                    storeOverflow(m,k:end) = 0;
                end
                break
            end
        else

        %% Limit Cycle Oscillation
        thetaStep = pi/2;
        thetaDotStep = 0;
        for i = 1:nStep
        
            %% Gaining Kineteic Energy
            if stage == 1 % Gaining kinetic energy need to Squat
        
                % Ending System Properties
                endCond = [0 1];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Squat; Squat],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                % Adjust Values
                timeDamp = timeN(1:end)+currentTime;
                currentTime = timeDamp(end);
                timeN(end+1) = timeN(end);
                thetaDampN(end+1) = thetaDampN(end);
                [thetaDotInstant] = conserveMomentum(thetaDotDampN(end),Swing,Squat,Stand);
                thetaDotDampN(end+1) = thetaDotInstant;

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:end),thetaDotDampN(1:end),'g')
                end

                % Store Values
                thetaStep(i+1,:) = real(thetaDampN(end));
                thetaDotStep(i+1,:) = real(thetaDotDampN(end));
                storeState(end+1:end+size(timeDamp',1)-1,1) = 0; % Stationary
                storeTime(end+1:end+size(timeDamp',1)-1,1) = timeDamp(2:end)';
                storeTheta(end+1:end+length(thetaDampN)-1,1) = thetaDampN(2:end);
                storeThetaDot(end+1:end+length(thetaDotDampN)-1,1) = thetaDotDampN(2:end);
        
                nrgTotal(i+1) = nrgT(thetaDampN(end),thetaDotDampN(end),Ist,Lst,Hst);
                nrgCheck = nrgTotal(i+1) - nrgT(0,0,Ist,Lst,Hst);
                stage = 2;
        
            %% Gaining Potential Energy
            elseif stage == 2 % Gaining potential energy need to Stand
                
                % Ending System Properties
                endCond = [0 2];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Stand; Stand],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                % Adjust time
                timeDamp = timeN(1:end)+currentTime;
                currentTime = timeDamp(end);

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:end),thetaDotDampN(1:end),'g')
                end
        
                % Store Values
                thetaStep(i+1,:) = real(thetaDampN(end));
                thetaDotStep(i+1,:) = real(thetaDotDampN(end));
                storeState(end+1:end+size(timeDamp',1)-1,1) = 0; % Stationary
                storeTime(end+1:end+size(timeDamp',1)-1,1) = timeDamp(2:end)';
                storeTheta(end+1:end+length(thetaDampN)-1,1) = thetaDampN(2:end);
                storeThetaDot(end+1:end+length(thetaDotDampN)-1,1) = thetaDotDampN(2:end);
        
                nrgTotal(i+1) = nrgT(thetaDampN(end),thetaDotDampN(end),Isq,Lsq,Hsq);
                nrgCheck = nrgTotal(i+1)-nrgT(0,0,Isq,Lsq,Hsq);
                stage = 1;
            end
        end

        disp('Limit Cycle Oscillation')
        overflow = max(abs(thetaDampN));
        end
        end

        storeOverflow(l,k) = overflow;
        end
        clf
    end
end

overInd = find(storeOverflow>pi);
storeOverflow(overInd) = pi;

figure
surf(Damping,Mass,storeOverflow)
shading interp

% a = gca();                      
% a.YScale = 'log';               % Set Zscale to logorithmic (comment out for linear)

xlabel('Damping');
ylabel('Mass');
zlabel('Theta')

toc
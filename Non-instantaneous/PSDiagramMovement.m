%% Kiiking Optimized Diagrams
% Cody Langston
clc; clear all; close all;
addpath(fullfile(pwd, 'NumericalProcesses'));
addpath(fullfile(pwd, 'Figures'));
pool = parpool("Threads");

Mass = linspace(5,0,20);
Length = linspace(0.1,5,10);

Damping = linspace(0,0.05,30); % Viscous Damping
dampTheta = pi/6;             % Viscous Damping
% Damping = linspace(0,2,40);    % Coulomb Damping
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
        Ms = 0.157;      % Mass of swing
        g = 9.8;         % Gravity
        
        Mp = Mass(l);   % Mass of person
        Ls = 0.406;     % Length of swing
        Lsq = 0.305;    % Length of squat
        Lst = 0.29;     % Length of stand

        % Mp = 0.3;
        % Ls = Length(l);
        % Lsq = Ls;
        % Lst = Ls-0.015;
        
        datumPE = Ls; % PE Datum
        Hs = datumPE-Ls/2;
        Hsq = datumPE-Lsq;
        Hst = datumPE-Lst;
        
        % Moment of inertia
        Is = 1/3*Ms*Ls^2;   % Moment of inertia swing
        Isq = Mp*Lsq^2;     % Moment of inertia squat
        Ist = Mp*Lst^2;     % Moment of inertia stand

        omegaN = sqrt(g*(Ms*Ls/2+Mp*Lsq)/(Is+Isq));
        mu = 2*Damping(k)*omegaN;
        
        % Initial Conditions
        thetaI = pi/2;              % Initial angle (rad)
        thetaDotI = 0;              % Initial velocity (rad/s)
        
        % Time parameters
        odeTime = 0.8;          % Simulation and overflow time for ODE45
        nStep = 70;                % Number of total cycles
        n = 20;                    % Number of undamped values per cycle
        moveTime = 0.2;
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
        for i = 1:9
        
            %% Gaining Kineteic Energy
            if stage == 1 % Gaining kinetic energy need to squat
        
                % Ending System Properties
                endCond = [0 1];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Squat; Squat],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                parfor j = 1:n
                    [tempTheta,tempThetaDot,~] = runODE(thetaDampN(j),thetaDotDampN(j),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Damped Ideal
        
                    nrgMove(j) = nrgT(tempTheta(end),tempThetaDot(end),Ist,Lst,Hst)
                end
        
                [~, indMove] = max(nrgMove);
                [thetaMove,thetaDotMove,timeMove] = runODE(thetaDampN(indMove),thetaDotDampN(indMove),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Damped Ideal
        
                % Adjust time
                timeDamp = timeN(1:indMove)+currentTime;
                currentTime = timeDamp(end);
                timeMove = timeMove+currentTime;
                currentTime = timeMove(end);

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:indMove),thetaDotDampN(1:indMove),'g')
                    plot(thetaMove,thetaDotMove,'r');
                end

                % Store Values
                thetaStep(i+1,:) = real(thetaMove(end));
                thetaDotStep(i+1,:) = real(thetaDotMove(end));
                storeState(end+1:end+size(timeDamp',1)-1,1) = 0; % Stationary
                storeState(end+1:end+size(timeMove,1)-1,1) = 2; % Squatting
                storeTime(end+1:end+size(timeDamp',1)-1,1) = timeDamp(2:end)';
                storeTime(end+1:end+size(timeMove,1)-1,1) = timeMove(2:end);
                storeTheta(end+1:end+indMove-1,1) = thetaDampN(2:indMove);
                storeTheta(end+1:end+size(thetaMove,1)-1,1) = thetaMove(2:end);
                storeThetaDot(end+1:end+indMove-1,1) = thetaDotDampN(2:indMove);
                storeThetaDot(end+1:end+size(thetaDotMove,1)-1,1) = thetaDotMove(2:end);

                if thetaDotMove(end)*thetaMove(end)<0
                    [thetaZero,thetaDotZero,timeZero] = runODE(thetaMove(end),thetaDotMove(end),moveInput,Swing,[Stand; Stand],odeTime,endCond);
                    if plotFigure
                        plot(thetaZero,thetaDotZero,'g')
                    end

                    thetaStep(i+1,:) = real(thetaZero(end));
                    thetaDotStep(i+1,:) = real(thetaDotZero(end));
                    storeState(end+1:end+size(timeZero,1)-1,1) = 0; % Stationary
                    storeTime(end+1:end+size(timeZero,1)-1,1) = timeZero(2:end);
                    storeTheta(end+1:end+size(thetaZero,1)-1,1) = thetaZero(2:end);
                    storeThetaDot(end+1:end+size(thetaDotZero,1)-1,1) = thetaDotZero(2:end);
                end

                nrgTotal(i+1) = nrgT(thetaDampN(indMove),thetaDotDampN(indMove),Ist,Lst,Hst);
                nrgCheck = nrgTotal(i+1) - nrgT(0,0,Ist,Lst,Hst);
                stage = 2;
        
            %% Gaining Potential Energy
            elseif stage == 2 % Gaining potential energy need to stand
                
                % Ending System Properties
                endCond = [0 2];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Stand; Stand],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                parfor j = 1:n
                    [tempTheta,tempThetaDot,~] = runODE(thetaDampN(j),thetaDotDampN(j),moveInput,Swing,[Stand; Squat],moveTime,[0 3]);   % Damped Ideal
        
                    nrgMove(j) = nrgT(tempTheta(end),tempThetaDot(end),Ist,Lst,Hst)
                end
        
                [~, indMove] = max(nrgMove);
                [thetaMove,thetaDotMove,timeMove] = runODE(thetaDampN(indMove),thetaDotDampN(indMove),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Damped Ideal
        
                % Adjust time
                timeDamp = timeN(1:indMove)+currentTime;
                currentTime = timeDamp(end);
                timeMove = timeMove+currentTime;
                currentTime = timeMove(end);

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:indMove),thetaDotDampN(1:indMove),'g')
                    plot(thetaMove,thetaDotMove,'b');
                end
        
                % Store Values
                timeDamp = timeDamp';
                thetaStep(i+1,:) = real(thetaMove(end));
                thetaDotStep(i+1,:) = real(thetaDotMove(end));
                storeState(end+1:end+size(timeDamp,1)-1,1) = 0; % Stationary
                storeState(end+1:end+size(timeMove,1)-1,1) = 1; % Standing
                storeTime(end+1:end+size(timeDamp,1)-1,1) = timeDamp(2:end);
                storeTime(end+1:end+size(timeMove,1)-1,1) = timeMove(2:end);
                storeTheta(end+1:end+indMove-1,1) = thetaDampN(2:indMove);
                storeTheta(end+1:end+size(thetaMove,1)-1,1) = thetaMove(2:end);
                storeThetaDot(end+1:end+indMove-1,1) = thetaDotDampN(2:indMove);
                storeThetaDot(end+1:end+size(thetaDotMove,1)-1,1) = thetaDotMove(2:end);

                if thetaDotMove(end)*thetaMove(end)>0
                    [thetaZero,thetaDotZero,timeZero] = runODE(thetaMove(end),thetaDotMove(end),moveInput,Swing,[Stand; Stand],odeTime,endCond);
                    if plotFigure
                        plot(thetaZero,thetaDotZero,'g')
                    end

                    thetaStep(i+1,:) = real(thetaZero(end));
                    thetaDotStep(i+1,:) = real(thetaDotZero(end));
                    storeState(end+1:end+size(timeZero,1)-1,1) = 0; % Stationary
                    storeTime(end+1:end+size(timeZero,1)-1,1) = timeZero(2:end);
                    storeTheta(end+1:end+size(thetaZero,1)-1,1) = thetaZero(2:end);
                    storeThetaDot(end+1:end+size(thetaDotZero,1)-1,1) = thetaDotZero(2:end);
                end
        
                nrgTotal(i+1) = nrgT(thetaDampN(indMove),thetaDotDampN(indMove),Isq,Lsq,Hsq);
                nrgCheck = nrgTotal(i+1)-nrgT(0,0,Isq,Lsq,Hsq);
                stage = 1;
            end
        end
            
        if nrgTotal(9) >= nrgTotal(7)
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
            if stage == 1 % Gaining kinetic energy need to squat
        
                % Ending System Properties
                endCond = [0 1];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Squat; Squat],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                parfor j = 1:n
                    [tempTheta,tempThetaDot,~] = runODE(thetaDampN(j),thetaDotDampN(j),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Damped Ideal
        
                    nrgMove(j) = nrgT(tempTheta(end),tempThetaDot(end),Ist,Lst,Hst)
                end
        
                [~, indMove] = max(nrgMove);
                [thetaMove,thetaDotMove,timeMove] = runODE(thetaDampN(indMove),thetaDotDampN(indMove),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Damped Ideal
        
                % Adjust time
                timeDamp = timeN(1:indMove)+currentTime;
                currentTime = timeDamp(end);
                timeMove = timeMove+currentTime;
                currentTime = timeMove(end);

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:indMove),thetaDotDampN(1:indMove),'g')
                    plot(thetaMove,thetaDotMove,'r');
                end

                % Store Values
                thetaStep(i+1,:) = real(thetaMove(end));
                thetaDotStep(i+1,:) = real(thetaDotMove(end));
                storeState(end+1:end+size(timeDamp',1)-1,1) = 0; % Stationary
                storeState(end+1:end+size(timeMove,1)-1,1) = 2; % Squatting
                storeTime(end+1:end+size(timeDamp',1)-1,1) = timeDamp(2:end)';
                storeTime(end+1:end+size(timeMove,1)-1,1) = timeMove(2:end);
                storeTheta(end+1:end+indMove-1,1) = thetaDampN(2:indMove);
                storeTheta(end+1:end+size(thetaMove,1)-1,1) = thetaMove(2:end);
                storeThetaDot(end+1:end+indMove-1,1) = thetaDotDampN(2:indMove);
                storeThetaDot(end+1:end+size(thetaDotMove,1)-1,1) = thetaDotMove(2:end);

                if thetaDotMove(end)*thetaMove(end)<0
                    [thetaZero,thetaDotZero,timeZero] = runODE(thetaMove(end),thetaDotMove(end),moveInput,Swing,[Stand; Stand],odeTime,endCond);
                    if plotFigure
                        plot(thetaZero,thetaDotZero,'g')
                    end

                    thetaStep(i+1,:) = real(thetaZero(end));
                    thetaDotStep(i+1,:) = real(thetaDotZero(end));
                    storeState(end+1:end+size(timeZero,1)-1,1) = 0; % Stationary
                    storeTime(end+1:end+size(timeZero,1)-1,1) = timeZero(2:end);
                    storeTheta(end+1:end+size(thetaZero,1)-1,1) = thetaZero(2:end);
                    storeThetaDot(end+1:end+size(thetaDotZero,1)-1,1) = thetaDotZero(2:end);
                end
        
                nrgTotal(i+1) = nrgT(thetaDampN(indMove),thetaDotDampN(indMove),Ist,Lst,Hst);
                nrgCheck = nrgTotal(i+1) - nrgT(0,0,Ist,Lst,Hst);
                stage = 2;
        
            %% Gaining Potential Energy
            elseif stage == 2 % Gaining potential energy need to stand
                
                % Ending System Properties
                endCond = [0 2];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Stand; Stand],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                parfor j = 1:n
                    [tempTheta,tempThetaDot,~] = runODE(thetaDampN(j),thetaDotDampN(j),moveInput,Swing,[Stand; Squat],moveTime,[0 3]);   % Damped Ideal
        
                    nrgMove(j) = nrgT(tempTheta(end),tempThetaDot(end),Ist,Lst,Hst)
                end
        
                [~, indMove] = max(nrgMove);
                [thetaMove,thetaDotMove,timeMove] = runODE(thetaDampN(indMove),thetaDotDampN(indMove),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Damped Ideal
        
                % Adjust time
                timeDamp = timeN(1:indMove)+currentTime;
                currentTime = timeDamp(end);
                timeMove = timeMove+currentTime;
                currentTime = timeMove(end);

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:indMove),thetaDotDampN(1:indMove),'g')
                    plot(thetaMove,thetaDotMove,'b');
                end
        
                % Store Values
                timeDamp = timeDamp';
                thetaStep(i+1,:) = real(thetaMove(end));
                thetaDotStep(i+1,:) = real(thetaDotMove(end));
                storeState(end+1:end+size(timeDamp,1)-1,1) = 0; % Stationary
                storeState(end+1:end+size(timeMove,1)-1,1) = 1; % Standing
                storeTime(end+1:end+size(timeDamp,1)-1,1) = timeDamp(2:end);
                storeTime(end+1:end+size(timeMove,1)-1,1) = timeMove(2:end);
                storeTheta(end+1:end+indMove-1,1) = thetaDampN(2:indMove);
                storeTheta(end+1:end+size(thetaMove,1)-1,1) = thetaMove(2:end);
                storeThetaDot(end+1:end+indMove-1,1) = thetaDotDampN(2:indMove);
                storeThetaDot(end+1:end+size(thetaDotMove,1)-1,1) = thetaDotMove(2:end);

                if thetaDotMove(end)*thetaMove(end)>0
                    [thetaZero,thetaDotZero,timeZero] = runODE(thetaMove(end),thetaDotMove(end),moveInput,Swing,[Stand; Stand],odeTime,endCond);
                    if plotFigure
                        plot(thetaZero,thetaDotZero,'g')
                    end

                    thetaStep(i+1,:) = real(thetaZero(end));
                    thetaDotStep(i+1,:) = real(thetaDotZero(end));
                    storeState(end+1:end+size(timeZero,1)-1,1) = 0; % Stationary
                    storeTime(end+1:end+size(timeZero,1)-1,1) = timeZero(2:end);
                    storeTheta(end+1:end+size(thetaZero,1)-1,1) = thetaZero(2:end);
                    storeThetaDot(end+1:end+size(thetaDotZero,1)-1,1) = thetaDotZero(2:end);
                end
        
                nrgTotal(i+1) = nrgT(thetaDampN(indMove),thetaDotDampN(indMove),Isq,Lsq,Hsq);
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
            if stage == 1 % Gaining kinetic energy need to squat
        
                % Ending System Properties
                endCond = [0 1];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Squat; Squat],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                parfor j = 1:n
                    [tempTheta,tempThetaDot,~] = runODE(thetaDampN(j),thetaDotDampN(j),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Damped Ideal
        
                    nrgMove(j) = nrgT(tempTheta(end),tempThetaDot(end),Ist,Lst,Hst)
                end
        
                [~, indMove] = max(nrgMove);
                [thetaMove,thetaDotMove,timeMove] = runODE(thetaDampN(indMove),thetaDotDampN(indMove),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Damped Ideal
        
                % Adjust time
                timeDamp = timeN(1:indMove)+currentTime;
                currentTime = timeDamp(end);
                timeMove = timeMove+currentTime;
                currentTime = timeMove(end);

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:indMove),thetaDotDampN(1:indMove),'g')
                    plot(thetaMove,thetaDotMove,'r');
                end

                % Store Values
                thetaStep(i+1,:) = real(thetaMove(end));
                thetaDotStep(i+1,:) = real(thetaDotMove(end));
                storeState(end+1:end+size(timeDamp',1)-1,1) = 0; % Stationary
                storeState(end+1:end+size(timeMove,1)-1,1) = 2; % Squatting
                storeTime(end+1:end+size(timeDamp',1)-1,1) = timeDamp(2:end)';
                storeTime(end+1:end+size(timeMove,1)-1,1) = timeMove(2:end);
                storeTheta(end+1:end+indMove-1,1) = thetaDampN(2:indMove);
                storeTheta(end+1:end+size(thetaMove,1)-1,1) = thetaMove(2:end);
                storeThetaDot(end+1:end+indMove-1,1) = thetaDotDampN(2:indMove);
                storeThetaDot(end+1:end+size(thetaDotMove,1)-1,1) = thetaDotMove(2:end);

                if thetaDotMove(end)*thetaMove(end)<0
                    [thetaZero,thetaDotZero,timeZero] = runODE(thetaMove(end),thetaDotMove(end),moveInput,Swing,[Stand; Stand],odeTime,endCond);
                    if plotFigure
                        plot(thetaZero,thetaDotZero,'g')
                    end

                    thetaStep(i+1,:) = real(thetaZero(end));
                    thetaDotStep(i+1,:) = real(thetaDotZero(end));
                    storeState(end+1:end+size(timeZero,1)-1,1) = 0; % Stationary
                    storeTime(end+1:end+size(timeZero,1)-1,1) = timeZero(2:end);
                    storeTheta(end+1:end+size(thetaZero,1)-1,1) = thetaZero(2:end);
                    storeThetaDot(end+1:end+size(thetaDotZero,1)-1,1) = thetaDotZero(2:end);
                end
        
                nrgTotal(i+1) = nrgT(thetaDampN(indMove),thetaDotDampN(indMove),Ist,Lst,Hst);
                nrgCheck = nrgTotal(i+1) - nrgT(0,0,Ist,Lst,Hst);
                stage = 2;
        
            %% Gaining Potential Energy
            elseif stage == 2 % Gaining potential energy need to stand
                
                % Ending System Properties
                endCond = [0 2];
                moveInput = [mu g moveTime];
                [thetaDamp,thetaDotDamp,timeDamp] = runODE(thetaStep(i,1),thetaDotStep(i,1),moveInput,Swing,[Stand; Stand],odeTime,endCond);   % Damped Ideal
        
                timeN = linspace(timeDamp(1),timeDamp(end),n);
                thetaDampN = interp1(timeDamp,thetaDamp,timeN);
                thetaDotDampN = interp1(timeDamp,thetaDotDamp,timeN);
        
                parfor j = 1:n
                    [tempTheta,tempThetaDot,~] = runODE(thetaDampN(j),thetaDotDampN(j),moveInput,Swing,[Stand; Squat],moveTime,[0 3]);   % Damped Ideal
        
                    nrgMove(j) = nrgT(tempTheta(end),tempThetaDot(end),Ist,Lst,Hst)
                end
        
                [~, indMove] = max(nrgMove);
                [thetaMove,thetaDotMove,timeMove] = runODE(thetaDampN(indMove),thetaDotDampN(indMove),moveInput,Swing,[Squat; Stand],moveTime,[0 3]);   % Damped Ideal
        
                % Adjust time
                timeDamp = timeN(1:indMove)+currentTime;
                currentTime = timeDamp(end);
                timeMove = timeMove+currentTime;
                currentTime = timeMove(end);

                % Phase Diagram
                if plotFigure
                    figure(1)
                    hold on
                    plot(thetaDamp,thetaDotDamp,'--k');
                    plot(thetaDampN(1:indMove),thetaDotDampN(1:indMove),'g')
                    plot(thetaMove,thetaDotMove,'b');
                end
        
                % Store Values
                timeDamp = timeDamp';
                thetaStep(i+1,:) = real(thetaMove(end));
                thetaDotStep(i+1,:) = real(thetaDotMove(end));
                storeState(end+1:end+size(timeDamp,1)-1,1) = 0; % Stationary
                storeState(end+1:end+size(timeMove,1)-1,1) = 1; % Standing
                storeTime(end+1:end+size(timeDamp,1)-1,1) = timeDamp(2:end);
                storeTime(end+1:end+size(timeMove,1)-1,1) = timeMove(2:end);
                storeTheta(end+1:end+indMove-1,1) = thetaDampN(2:indMove);
                storeTheta(end+1:end+size(thetaMove,1)-1,1) = thetaMove(2:end);
                storeThetaDot(end+1:end+indMove-1,1) = thetaDotDampN(2:indMove);
                storeThetaDot(end+1:end+size(thetaDotMove,1)-1,1) = thetaDotMove(2:end);

                if thetaDotMove(end)*thetaMove(end)>0
                    [thetaZero,thetaDotZero,timeZero] = runODE(thetaMove(end),thetaDotMove(end),moveInput,Swing,[Stand; Stand],odeTime,endCond);
                    if plotFigure
                        plot(thetaZero,thetaDotZero,'g')
                    end

                    thetaStep(i+1,:) = real(thetaZero(end));
                    thetaDotStep(i+1,:) = real(thetaDotZero(end));
                    storeState(end+1:end+size(timeZero,1)-1,1) = 0; % Stationary
                    storeTime(end+1:end+size(timeZero,1)-1,1) = timeZero(2:end);
                    storeTheta(end+1:end+size(thetaZero,1)-1,1) = thetaZero(2:end);
                    storeThetaDot(end+1:end+size(thetaDotZero,1)-1,1) = thetaDotZero(2:end);
                end
        
                nrgTotal(i+1) = nrgT(thetaDampN(indMove),thetaDotDampN(indMove),Isq,Lsq,Hsq);
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
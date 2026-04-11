function [thetaStore, thetaDotStore, tauStore, stateStore, thetaOptim] = runCycle(thetaI,state)

%% Initialize
n = 100;
odeTime = 0.75*2*pi;
thetaDotI = 0;
g = 9.8;

sigma = state(1);   % Nondimensional Mass Mr/Ms
delta = state(2);   % Nondimensional Length 1-Lst/Ls
beta = state(3);    % Nondimensional Time Tst/Tperiod*2*pi
zeta = state(4);    % Nondimensional Damping

% Energy Calculations
nrgTsq = @(theta,thetaDot) 1/6*thetaDot^2+1/3*(1-cos(theta));

%% Squat
endCond = [0 1];        % Stop at theta = 0
moveInput = 1;          % Use squat equations
state = [sigma 0 odeTime zeta];
[thetaTemp,thetaDotTemp,tauTemp] = runNDODE(thetaI,thetaDotI,state,moveInput,endCond);

tauN = linspace(tauTemp(1),tauTemp(end),n);
thetaN = interp1(tauTemp,thetaTemp,tauN);
thetaDotN = interp1(tauTemp,thetaDotTemp,tauN);

% plot(thetaDampN,thetaDotDampN);
% hold on

%% Standing
endCond = [0 2];        % Stop at time limit
moveInput = 1;          % Use standing equations
state = [sigma delta beta zeta];

% Optimize standing time
nA = 1;     % lowerBound
nD = n-1;   % upperBound
nrgMoveA = 0;
nrgMoveD = 0;
while nD-nA >= 4
    nB = floor(nA+0.33*(nD-nA));    % boundC = nA + 1/3(nD-nA)
    nC = floor(nA+0.66*(nD-nA));    % boundD = nA + 2/3(nD-nA)
    
    [thetaB,thetaDotB,tauB] = runNDODE(thetaN(nB),thetaDotN(nB),state,moveInput,endCond);
    nrgMoveB = nrgTsq(thetaB(end),thetaDotB(end));
    [thetaC,thetaDotC,tauC] = runNDODE(thetaN(nC),thetaDotN(nC),state,moveInput,endCond);
    nrgMoveC = nrgTsq(thetaC(end),thetaDotC(end));

    % plot([nA, nB, nC, nD], [nrgMoveA, nrgMoveB, nrgMoveC, nrgMoveD])

    if nrgMoveB > nrgMoveC      % If maximum between nA and nC
        nD = nC;
        nrgMoveD = nrgMoveC;
    else                        % If maximum between nB and ND
        nA = nB;
        nrgMoveA = nrgMoveB;
        nMax = nC;
    end
end

% Store optimal time
if nrgMoveB > nrgMoveC      % If maximum between nA and nC
    nMax = nC;
    thetaOptim = thetaN(nMax);
    thetaStore = [thetaN(1:nMax-1), thetaC'];
    thetaDotStore = [thetaDotN(1:nMax-1), thetaDotC'];
    tauStore = [tauN(1:nMax-1), tauC'+tauN(nMax)];
    stateStore = [ones(1,nMax-1), 2*ones(size(tauC'))];
else                        % If maximum between nB and ND
    nMax = nB;
    thetaOptim = thetaN(nMax);
    thetaStore = [thetaN(1:nMax-1), thetaB'];
    thetaDotStore = [thetaDotN(1:nMax-1), thetaDotB'];
    tauStore = [tauN(1:nMax-1), tauB'+tauN(nMax)];
    stateStore = [ones(1,nMax-1), 2*ones(size(tauB'))];
end

%% Stand
endCond = [0 2];        % Stop at thetaDot = 0
moveInput = 2;          % Use stand equations
state = [sigma 0 odeTime zeta];
[thetaTemp,thetaDotTemp,tauTemp] = runNDODE(thetaStore(end),thetaDotStore(end),state,moveInput,endCond);

thetaStore = [thetaStore(1:end-1), thetaTemp'];
thetaDotStore = [thetaDotStore(1:end-1), thetaDotTemp'];
tauStore = [tauStore(1:end-1), tauTemp'+tauStore(end)];
stateStore = [stateStore(1:end-1), 3*ones(1,size(tauTemp,1))];







end




% parfor j = 1:n
%     [thetaTemp,thetaDotTemp,~] = runNDODE(thetaDampN(j),thetaDotDampN(j),state,moveInput,endCond);   % Damped Ideal
%     nrgMove(j) = nrgTsq(thetaTemp(end),thetaDotTemp(end),sigma);
% end
% figure
% plot(nrgMove)
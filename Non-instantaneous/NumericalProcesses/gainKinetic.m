function [theta,thetaDot] = gainKinetic(thetaI,thetaDotI,thetaF,g,mu,swing,state,n)

    %% Store Variables
    Ms = swing(1);
    Ls = swing(2);
    Hs = swing(3);
    Is = swing(4);

    Mp = state(1);
    Lp = state(2);
    Hp = state(3);
    Ip = state(4);

    %% Energy
    nrgKE = @(thetaDot,I) 1/2*I*thetaDot.^2;               % Rotational Kinetic Energy
    nrgPE = @(theta,M,L,H) M*g*(H+L*(1-cos(theta)));        % Gravitational Potential Energy
    % Total Energy
    nrgT = @(theta,thetaDot) (nrgPE(theta,Mp,Lp,Hp)+nrgPE(theta,Ms,Ls/2,Hs))...  % Potential
         + (nrgKE(thetaDot,Ip)+nrgKE(thetaDot,Is));                                   % Rotational Kinetic

    %% Calcualtions
    theta = linspace(thetaI,thetaF,n);

    % Initial Energy
    nrgTI = nrgT(thetaI,thetaDotI);

    % Final Potential Energy
    nrgPEf = nrgPE(theta,Mp,Lp,Hp)+nrgPE(theta,Ms,Ls/2,Hs);
    nrgKEf = nrgTI-nrgPEf;

    % Final Angular Velocity
    thetaDot = sqrt((nrgKEf)/(1/2*(Ip+Is)));

    %% Damping Calculations
    nrgDamp = cumtrapz(theta,-mu*thetaDot);

    nrgKEf = nrgTI-nrgPEf-abs(nrgDamp);

    thetaDot = sqrt((nrgKEf)/(1/2*(Ip+Is)));

    % Direction
    if thetaI>=0
        thetaDot = -thetaDot;
    end
end
function [theta,thetaDot] = gainPotential(thetaI,thetaDotI,thetaDotF,g,swing,state,n)

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
    thetaDot = linspace(thetaDotI,thetaDotF,n);

    % Initial Energy
    nrgTI = nrgT(thetaI,thetaDotI);

    % Final Kinetic Energy
    nrgKEf = nrgKE(thetaDot,Ip)+nrgKE(thetaDot,Is);
    nrgPEf = nrgTI-nrgKEf;

    % Final Angular Orientation
    theta = acos(1-((nrgPEf)-(Ms*g*Hs+Mp*g*Hp))/(Ms*g*Ls/2+Mp*g*Lp));

    % Direction
    if thetaDotI<=0
        theta = -theta;
    end
end
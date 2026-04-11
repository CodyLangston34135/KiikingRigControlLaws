function [thetaDot] = conserveMomentum(thetaDotI,swing,squat,stand)

    %% Store Variables
    Ms = swing(1);
    Ls = swing(2);
    Hs = swing(3);
    Is = swing(4);

    Msq = squat(1);
    Lsq = squat(2);
    Hsq = squat(3);
    Isq = squat(4);

    Mst = stand(1);
    Lst = stand(2);
    Hst = stand(3);
    Ist = stand(4);

    %% Calculations
    thetaDot = thetaDotI*(Isq+Is)/(Ist+Is);

end
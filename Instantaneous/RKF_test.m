%% Runge-Kutta-Fehlberg adaptive step size method 
tic
clear; close all;

% User Inputs
    epsilon = 0.001;     % accuracy threshold

% Initializing
    time = 16;
    t(1) = 0;
    x(1) = 0.01;
    i = 1;
    h = 0.0001;

    % define ODE
        f = @(t,x) x*(1 - x);

% Solution

c1 = 1/4; 
c2 = 3/8; c3 = 3/32; c4 = 9/32; 
c5 = 12/13; c6 = 1932/2197; c7 = -7200/2197; c8 = 7296/2197; 
c9 = 439/216; c10 = -8; c11 = 3680/513; c12 = -845/4104; 
c13 = 1/2; c14 = -8/27; c15 = 2; c16 = -3544/2565; c17 = 1859/4104; c18 = -11/40; 
c19 = 16/135; c20 = 6656/12825; c21 = 28561/56430; c22 = -9/50; c23 = 2/55;
c24 = 25/216; c25 = 1408/2565; c26 = 2197/4104; c27 = -1/5;     

    while t(i) < time
        K1 = h*f(t(i),x(i));
        K2 = h*f(t(i) + c1*h, x(i) + c1*K1);
        K3 = h*f(t(i) + c2*h, x(i) + c3*K1 + c4*K2);
        K4 = h*f(t(i) + c5*h, x(i) + c6*K1 + c7*K2 + c8*K3);
        K5 = h*f(t(i) + h, x(i) + c9*K1 + c10*K2 + c11*K3 + c12*K4);
        K6 = h*f(t(i) + c13*h, x(i) + c14*K1 + c15*K2 + c16*K3 + c17*K4 + c18*K5);
        RK5 = x(i) + c19*K1 + c20*K3 + c21*K4 + c22*K5 + c23*K6;
        RK4 = x(i) + c24*K1 + c25*K3 + c26*K4 + c27*K5;
        if abs(RK4 - RK5) < epsilon
            x(i+1) = RK5;
            t(i+1) = t(i) + h;
            h = 2*h;
            i = i + 1;
        else
            h = 0.5*h;
        end
    end

% define true solution
    T = linspace(0,16,1000);
    X = 1 - 1./(1 + (exp(T)./99));

    figure
    plot(T,X)
    hold on
    plot(t,x)
toc
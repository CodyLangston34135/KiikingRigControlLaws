%% Runge-Kutta-Fehlberg adaptive step size method 

function RKF(epsilon)
    for j = 1:length(epsilon)
        % Initializing
            time = 16;
            t(1) = 0;
            x(1) = 0.01;
            i = 1;
            h = 1;
        
            % define ODE
                f = @(t,x) x*(1 - x);
        
        % Solution 
        
            while t(i) < time
                K1 = h*f(t(i),x(i));
                K2 = h*f(t(i) + (1/4)*h, x(i) + (1/4)*K1);
                K3 = h*f(t(i) + (3/8)*h, x(i) + (3/32)*K1 + (9/32)*K2);
                K4 = h*f(t(i) + (12/13)*h, x(i) + (1932/2197)*K1 - (7200/2197)*K2 + (7296/2197)*K3);
                K5 = h*f(t(i) + h, x(i) + (439/216)*K1 - 8*K2 + (3680/513)*K3 - (845/4104)*K4);
                K6 = h*f(t(i) + (1/2)*h, x(i) - (8/27)*K1 + 2*K2 - (3544/2565)*K3 + (1859/4104)*K4 - (11/40)*K5);
                RK5 = x(i) + (16/135)*K1 + (6656/12825)*K3 + (28561/56430)*K4 - (9/50)*K5 + (2/55)*K6;
                RK4 = x(i) + (25/216)*K1 + (1408/2565)*K3 + (2197/4104)*K4 - (1/5)*K5;
                if abs(RK4 - RK5) < epsilon(j)
                    x(i+1) = RK5; %#ok<AGROW> 
                    t(i+1) = t(i) + h; %#ok<AGROW> 
                    h = 2*h;
                    i = i + 1;
                else
                    h = 0.5*h;
                end
            end
        
        % define true solution
            T = linspace(0,16,1000);
            X = 1 - 1./(1 + (exp(T)./99));
        
            figure(j)
            plot(T,X)
            hold on
            plot(t,x, 'o-')
            grid on
            ylabel('$x(t)$')
            xlabel('$t$')
            legend('Analyical Soltuion', 'Numerical Solution')
    end
end
function [xN, yOut] = testFit(x,y,xN,plotData)

% y = detrend(y);                                                                 % Remove Linear Trend

yu = max(y);
yl = min(y);
yr = (yu-yl);                                                                   % Range of ‘y’
yz = y-yu+(yr/2);
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                            % Returns Approximate Zero-Crossing Indices Of Argument Vector
zt = x(zci(y));
per = 2*mean(diff(zt));                                                         % Estimate period
ym = mean(y);                                                                   % Estimate offset

fit = @(b,x)  b(1) .* exp(-b(2).*x) .* (sin(2*pi*x./b(3) + 2*pi/b(4))) + b(5);  % Objective Function to fit
fcn = @(b) norm(fit(b,x) - y);                                                  % Least-Squares cost function
[s,nmrs] = fminsearch(fcn, [yr; -10;  per;  -1;  ym]);                           % Minimise Least-Squares

xp = linspace(min(x),max(x), 500);

if plotData
    plot(x,y,'b', 'LineWidth',1.5)
    hold on
    plot(xp,fit(s,xp), '--r')
    hold off
    grid
    xlabel('Time')
    ylabel('Amplitude')
    legend('Original Data',  'Fitted Curve')
    text(0.3*max(xlim),0.7*min(ylim), sprintf('$y = %.3f\\cdot e^{%.0f\\cdot x}\\cdot sin(2\\pi\\cdot x\\cdot %.0f%.3f)$', [s(1:2); 1./s(3:4)]), 'Interpreter','latex')
end

yOut = fit(s,xN);

end
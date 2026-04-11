function testSmooth(t,y)

Fs=60;
Fn = Fs/2;                                                                      % Nyquist Frequency
L = size(t,1);
% y = detrend(y);                                                                 % Remove Linear Trend

fty = fft(y)/L;
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                                             % Frequency Vector
Iv = 1:length(Fv);                                                              % Index Vector

figure(1)
plot(Fv, abs(fty(Iv)))
grid

Wp =   5/Fn;                                                                    % Passband (Normalised)
Ws =  20/Fn;                                                                    % Stopband (Normalised)
Rp =  1;                                                                        % Passband Ripple
Rs = 25;                                                                        % Stopband Ripple
[n,Wn]  = buttord(Wp,Ws,Rp,Rs);                                                 % Filter Order
[b,a]   = butter(n,Wn);                                                         % Transfer Function Coefficients
[sos,g] = tf2sos(b,a);                                                          % Second-Order-Section For Stability

figure(2)
freqz(sos, 2048, Fs)

yf = filtfilt(sos,g,y);

yu = max(yf);
yl = min(yf);
yr = (yu-yl);                                                                   % Range of ‘y’
yz = yf-yu+(yr/2);
zt = t(yz(:) .* circshift(yz(:),[1 0]) <= 0);                                   % Find zero-crossings
per = 2*mean(diff(zt));                                                         % Estimate period
ym = mean(y);                                                                   % Estimate offset

fit = @(b,x)  b(1) .* exp(b(2).*x) .* (sin(2*pi*x./b(3) + 2*pi/b(4))) + b(5);   % Objective Function to fit
fcn = @(b) sum((fit(b,t) - yf).^2);                                             % Least-Squares cost function
s = fminsearch(fcn, [yr; -1;  per;  -1;  ym]);                                   % Minimise Least-Squares

xp = linspace(min(t),max(t));

figure(3)
scatter(t,y,'b')
hold on
plot(t,yf,'y', 'LineWidth',2)
plot(xp,fit(s,xp), 'r', 'LineWidth',1.5)
hold off
grid
title('Sample Data 2')
xlabel('Time')
ylabel('Amplitude')
legend('Original Data', 'Lowpass-Filtered Data', 'Fitted Curve')
end
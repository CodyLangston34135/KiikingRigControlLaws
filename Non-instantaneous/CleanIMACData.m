%% CleanIMACData
clc; clear all; close all;

addpath(fullfile(pwd, 'NumericalProcesses'));
addpath(fullfile(pwd, 'ExperimentalProcesses'));
addpath(fullfile(pwd, 'ExperimentalData'));
addpath(fullfile(pwd, 'IMACData'));
addpath(fullfile(pwd, 'Figures'));

load('DampedMove.mat')

index = 1:size(thetaN,2);

figure(1)
subplot(3,1,1)
plot(thetaN)

% removeEnd = 3400;
removeEnd = size(thetaN,2);
removeData = [995:1013 removeEnd:size(thetaN,2)];

time = timeN;
indexN = index;
indexN(removeData) = [];
distN(removeData) = [];
thetaDotDotSN(removeData) = [];
thetaDotSN(removeData) = [];
thetaN(removeData) = [];
time(removeData) = [];
timeN = timeN(1:removeEnd-1);

subplot(3,1,2)
plot(indexN,thetaN)

distN = interp1(time,distN,timeN,"spline");
thetaDotDotSN = interp1(time,thetaDotDotSN,timeN,"spline");
thetaDotSN = interp1(time,thetaDotSN,timeN,"spline");
thetaN = interp1(time,thetaN,timeN,"spline");

subplot(3,1,3)
plot(thetaN);

save Clean timeN thetaN thetaDotSN thetaDotDotSN distN moveIndex


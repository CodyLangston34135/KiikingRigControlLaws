%% CalculateSledTime
clc; clear all; close all;

addpath(fullfile(pwd, 'NumericalProcesses'));
addpath(fullfile(pwd, 'ExperimentalProcesses'));
addpath(fullfile(pwd, 'ExperimentalData'));
addpath(fullfile(pwd, 'IMACData/WhiteBackground/Data'));
addpath(fullfile(pwd, 'Figures'));

load('UndampedMove.mat')

barLength = 13.5;

%% Format Standing/Squatting times
index = 1;
storeRun = [1,1,moveIndex(1)];
for i = 1:size(moveIndex,2)-1
    if moveIndex(i) == moveIndex(i+1)
        storeRun(index,2) = i+1;
    else
        index=index+1;
        storeRun(index,:) = [i+1, i+1, moveIndex(i+1)];
    end
end

for i = 1:size(storeRun,1)
    if storeRun(i,3) == 0
        distBounds(i) = mean(distN(storeRun(i,1):storeRun(i,2)));
    end
end
moveTime = (storeRun(:,2)-storeRun(:,1))*(timeN(2)-timeN(1));
moveTime = moveTime(1:2:end);
mean(moveTime)

figure
for i = 1:size(moveIndex,2)-1
    hold on
    if moveIndex(i)==1
        plot([timeN(i), timeN(i+1)],barLength*[distN(i), distN(i+1)],'b');
    elseif moveIndex(i)==2
        plot([timeN(i), timeN(i+1)],barLength*[distN(i), distN(i+1)],'r');
    else
        plot([timeN(i), timeN(i+1)],barLength*[distN(i), distN(i+1)],'g');
    end
end
title('Sled Movement')
ylabel('Pixel Distance')
xlabel('Time (s)')
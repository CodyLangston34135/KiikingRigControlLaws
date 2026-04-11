%% IMAC Video Creator
clc; clear all; close all;

addpath(fullfile(pwd, 'NumericalProcesses'));
addpath(fullfile(pwd, 'ExperimentalProcesses'));
addpath(fullfile(pwd, 'ExperimentalData'));
addpath(fullfile(pwd, 'IMACData/WhiteBackground/Data'));
addpath(fullfile(pwd, 'Figures'));

%% Initialize Data
% Load Data
dataInput = 'LCOMove.mat';
load(dataInput);

% Load Video
dataDir = './Data';
filename = 'LCOMove.MP4';
vidFile = fullfile(dataDir, filename);
vr = VideoReader(vidFile);
[~, writeTag, ~] = fileparts(vidFile);
frameRate = vr.FrameRate; 
totalFrames = vr.NumFrames;

% Create Video Writer
v = VideoWriter('LCOMove.avi');
% v.FrameRate = 60;
open(v);

%% User Inputs
yBound = max(thetaN)+0.5;
barLength = 13.5;
periodStep = [22 45];

% Video Properties

% UndampedMove
% bounds = [1 2300];
% crop = [100 1050; 515, 1465]; % [yStart yEnd; xStart, xEnd]
% videoCut = [500 750 2444]; % [slowStart slowEnd vidEnd]
% vidSpeed = [4 1 4];

% DampedMove
% bounds = [1 1000];
% crop = [100 1050; 500, 1450]; % [yStart yEnd; xStart, xEnd]
% videoCut = [500 750 1250]; % [slowStart slowEnd vidEnd]
% vidSpeed = [4 1 4];

% Damped2MoveSpeedUp
% bounds = [3000 6850];
% crop = [100 1050; 550, 1500]; % [yStart yEnd; xStart, xEnd]
% videoCut = [4500 4750 6960]; % [slowStart slowEnd vidEnd]
% vidSpeed = [4 1 8];

% Damped2Move
% bounds = [3000 5868];
% crop = [100 1050; 550, 1500]; % [yStart yEnd; xStart, xEnd]
% videoCut = [4500 4750 5868]; % [slowStart slowEnd vidEnd]
% vidSpeed = [4 1 4];

% LCOMove
bounds = [1500 4368];
crop = [100 1050; 425, 1375]; % [yStart yEnd; xStart, xEnd]
videoCut = [3000 3250 200]; % [slowStart slowEnd vidEnd]
vidSpeed = [4 1 4];

%% Create Plots
% Format thetaN
thetaN = thetaN-mean(thetaN);

if bounds(1)~=1
    timeN(bounds(1):bounds(2)) = timeN(1:bounds(2)-bounds(1)+1);
end

slowInd1 = bounds(1):vidSpeed(1):videoCut(1);
slowInd2 = videoCut(1)+1:vidSpeed(2):videoCut(2);
slowInd3 = videoCut(2)+1:vidSpeed(3):bounds(2)-periodStep(1);

slowInd = [slowInd1, slowInd2, slowInd3];
slowPeriod = [slowInd1+periodStep(2), slowInd2+periodStep(2), slowInd3+periodStep(2)];

vidLength = size(slowInd,2)/v.FrameRate + size(bounds(2)-periodStep(1)+2:vidSpeed(3):videoCut(3),2)/v.FrameRate

% Plot Angular Displacement
figure
% set(gcf,'position',[0, 0, crop(2,2)-crop(2,1)+1, crop(1,2)-crop(1,1)+1])
set(gcf,'position',[0, 0, 1000, crop(1,2)-crop(1,1)+1])

subplot(2,1,1)
subplot(2,1,2)
g = gcf;
a1 = g.Children(1);
a1.Position = [(1-.45)/2 .1 .45 .45];

a2 = g.Children(2);
a2.Position = [(1-.8)/2 .65 .8 .3];

for i = 1:size(slowInd,2)-1

    % Get Video
    vid = read(vr, slowInd(i));
    vid = vid(crop(1,1):crop(1,2),crop(2,1):crop(2,2),:,:);

    % Create Figure
    axes(a2)
    ylabel('Angular Displacement (rad)')
    xlabel('Time (s)')
    xlim([0 timeN(bounds(2))])
    ylim([-yBound yBound])
    hold on
    box on
    plot(timeN(slowInd(i):slowInd(i+1)),thetaN(slowInd(i):slowInd(i+1)),'b','LineWidth',2);
    
    axes(a1)
    plot(thetaN(slowInd(i):slowInd(i+1)),thetaN(slowPeriod(i):slowPeriod(i+1)),'b','LineWidth',1.5);
    hold on
    ylabel('Angular Displacement (rad)')
    xlabel('Angular Displacement (rad)')
    xlim([-yBound yBound])
    ylim([-yBound yBound])


    frame = getframe(gcf);

    frame = uint8(frame.cdata);

    % Mesh Video and Figure
    blank = zeros([size(vid,1),size(vid,2)+size(frame,2),3]);
    blank = uint8(blank);
    blank(1:size(vid,1),1:size(vid,2),:) = vid;
    blank(1:size(vid,1),size(vid,2)+1:end,:) = frame;

    writeVideo(v,blank);
end

if videoCut(3) ~= bounds(2)
    for i = bounds(2)-periodStep(1)+2:vidSpeed(3):videoCut(3)

    % Get Video
    vid = read(vr, i);
    vid = vid(crop(1,1):crop(1,2),crop(2,1):crop(2,2),:,:);

    % Mesh Video and Figure
    blank = zeros([size(vid,1),size(vid,2)+size(frame,2),3]);
    blank = uint8(blank);
    blank(1:size(vid,1),1:size(vid,2),:) = vid;
    blank(1:size(vid,1),size(vid,2)+1:end,:) = frame;

    writeVideo(v,blank);
    end
end

close(v)


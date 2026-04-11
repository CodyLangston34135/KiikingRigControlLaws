%% Experimental Phase Diagrams
% Cody Langston
clc; clear all; close all;
tic

addpath(fullfile(pwd, 'Data'));
addpath(fullfile(pwd, 'ExperimentalProcesses'));
addpath(fullfile(pwd, 'Figures'));
% pool = parpool("Threads");

%% User Inputs
dataDir = './Data';
filename = 'LCOMove.MP4';

% Define datum
% datumPos = [505 442];                                           % Position of axle
datumPos = [510 430];
datumDistRed = [375 400];                                       % Estimated distance for red dot
datumDistBlue = 430;                                            % Estimated distance for blue dot
barLength = 13.5;                                               % Length of bar
datumDist = [datumDistRed datumDistBlue];

% Define allowed changes to account for noise
changeTheta = 900;
changeDistRed = 60;
changeDistBlue = 45;
maxChange = [changeTheta changeDistRed changeDistBlue];

% Define video processing setup
segments = 10;
periodSteps = 16;

showFigure = false;
storeVideo = false;
storeData = true;
debugTheta = false;

vidFile = fullfile(dataDir, filename);

% cd "G:\DDSL Lab\Johannes Kiiking\MovementProcessing\Data"
% ffmpeg -i input.mp4 -vcodec libx264 output.mp4


%% Load Video
vr = VideoReader(vidFile);
[~, writeTag, ~] = fileparts(vidFile);
FrameRate = vr.FrameRate; 
totalFrames = vr.NumFrames;
% totalFrames = 1000;
videoTime = floor(totalFrames/FrameRate);
fprintf('Video Length %i Seconds, %i Frames\n',videoTime,totalFrames)

% Segment frames to use
useFrames = floor(linspace(0,totalFrames,segments));

% Initialize loop variables
storeTheta = zeros(1,0);
storeBlueDist = zeros(1,0);
storeRedDist = zeros(1,0);
storeVid = zeros(1,1,1,0);

for j = 1:segments-1
    % Read and crop video
    vid = read(vr, [useFrames(j)+1,useFrames(j+1)]);
    vid = vid(100:1050,425:1375,:,:);         % Output15

    % Mask
    % vid(330:550,480:650,:,:) = 0;
    % vid(1:200,1:100,:,:) = 0;
    % vid(35:80,260:300,:,:) = 0;
    % vid(79:140,273:320,:,:) = 0;
    % vid(5:50,700:840,:,:) = 0;
    % vid(830:840,630:645,:,:) = 0;

    for i = 1:size(vid,4)                   % PARFOR
        % Store video frame and swap to YCbCr domain
        tempRGB = vid(:,:,:,i);
        tempYCbCr = rgb2ycbcr(tempRGB);
    
        % Find regions with lots of blue and red
        blue = tempYCbCr(:,:,2)>=150;
        blue = bwareaopen(blue,10);

        red = tempYCbCr(:,:,3)>=140;
        red = bwareaopen(red,10);

        % Store if necessary
        if storeVideo
            vidBlank(:,:,:,i) = blue+red;
        end

        % Show thresholds if necessary
        if showFigure
            subplot(1,3,1)
            imshow(tempRGB)
            impixelinfo
            subplot(1,3,2)
            imshow(blue)
            impixelinfo
            subplot(1,3,3)
            imshow(red)
            impixelinfo
        end
    
        %% Process Blue
        % Label blue regions
        [blueLabel, bluePoints] = bwlabel(blue);
        s = regionprops(blueLabel, 'Centroid');
    
        % If no blue regions dont store anything
        if bluePoints == 0
            theta(i) = 0;
            blueDist(i) = 0;
        else
            % Store centroid
            tempCent = [0 0];
            for k = 1:bluePoints
                tempCent(k,:) = s(k).Centroid;
            end

            % Calculate distance from datum
            tempDist = sqrt((tempCent(:,1)-datumPos(:,1)).^2+(tempCent(:,2)-datumPos(:,2)).^2); % Distance Formula
            [~, ind] = min(tempDist-datumDistBlue);
    
            % Store values
            blueCent = s(ind).Centroid;  % [x y]
            blueDist(i) = tempDist(ind);
            theta(i) = atan2d(-(blueCent(1)-datumPos(1)),(blueCent(2)-datumPos(2))); % Angle off horizontal
        end

        %% Process Red
        % Label red regions
        [redLabel, redPoints] = bwlabel(red);
        s = regionprops(redLabel, 'Centroid');
    
        % If no red regions dont store anything
        if redPoints == 0
            redDist(i) = 0;
        else
            % Store centrid
            tempCent = [0 0];
            for k = 1:redPoints
                tempCent(k,:) = s(k).Centroid;
            end

            % Calculate distance from datum
            tempDist = sqrt((tempCent(:,1)-datumPos(:,1)).^2+(tempCent(:,2)-datumPos(:,2)).^2); % Distance Formula
            [~, ind] = min(abs(tempDist-datumDistRed(1)));
    
            % Store values
            redDist(i) = tempDist(ind);
        end
    end
    storeTheta(1,end+1:end+size(theta,2)) = theta;
    storeBlueDist(1,end+1:end+size(blueDist,2)) = blueDist;
    storeRedDist(1,end+1:end+size(redDist,2)) = redDist;

    if storeVideo
        storeVid(1:size(vidBlank,1),1:size(vidBlank,2),:,end+1:end+size(theta,2)) = vidBlank;
    end
    fprintf('Finished Segment %i\n',j)
end

toc

%% Noise Reduction
indexN = 1:length(storeTheta);
timeN = indexN/FrameRate;
[theta,thetaDot,blueDist,redDist,time] = formatData(storeTheta,storeBlueDist,storeRedDist,timeN,datumDist,maxChange,debugTheta);

timeInd = find(timeN == time(end));
indexN = 1:timeInd;
timeN = indexN/FrameRate;

thetaDotS = smooth(thetaDot);
thetaDotS = thetaDotS';

thetaN = interp1(time,theta,timeN,"spline");
thetaDotSN = interp1(time,thetaDotS,timeN,"spline");
blueDistN = interp1(time,blueDist,timeN,"spline");
redDistN = interp1(time,redDist,timeN,"spline");

%% Sled movement calculations

% Calculate sled movement
distN = redDistN./blueDistN;
velN = gradient(distN,timeN);
[moveIndex,moveLength] = formatDist(distN,velN); % 0 = stationary, 1 = standing, 2 = squatting

moveTime = 1/FrameRate*mean(moveLength);

%% Create Plots

% Convert to radians
thetaN = thetaN*pi/180;
thetaDotSN = thetaDotSN*pi/180;

% Plot Angular Displacement
figure
subplot(3,1,1)
for i = 1:size(moveIndex,2)-1
    hold on
    if moveIndex(i)==1
        plot([timeN(i), timeN(i+1)],[thetaN(i), thetaN(i+1)],'b');
    elseif moveIndex(i)==2
        plot([timeN(i), timeN(i+1)],[thetaN(i), thetaN(i+1)],'r');
    else
        plot([timeN(i), timeN(i+1)],[thetaN(i), thetaN(i+1)],'g');
    end
end
title('Angular Displacement')
ylabel('Angular Displacement (rad)')
xlabel('Time (s)')
xlim([0 timeN(end)])

% Plot Angular Velocity
subplot(3,1,2)
for i = 1:size(moveIndex,2)-1
    hold on
    if moveIndex(i)==1
        plot([timeN(i), timeN(i+1)],[thetaDotSN(i), thetaDotSN(i+1)],'b');
    elseif moveIndex(i)==2
        plot([timeN(i), timeN(i+1)],[thetaDotSN(i), thetaDotSN(i+1)],'r');
    else
        plot([timeN(i), timeN(i+1)],[thetaDotSN(i), thetaDotSN(i+1)],'g');
    end
end
title('Angular Velocity')
ylabel('Angular Velocity (rad/s)')
xlabel('Time (s)')
xlim([0 timeN(end)])
ylim([-1500 1500]*pi/180)

% Plot Angular Acceleration
thetaDotDotN = gradient(thetaDotSN,timeN);
thetaDotDotSN = smooth(thetaDotDotN);
subplot(3,1,3)
for i = 1:size(moveIndex,2)-1
    hold on
    if moveIndex(i)==1
        plot([timeN(i), timeN(i+1)],[thetaDotDotSN(i), thetaDotDotSN(i+1)],'b');
    elseif moveIndex(i)==2
        plot([timeN(i), timeN(i+1)],[thetaDotDotSN(i), thetaDotDotSN(i+1)],'r');
    else
        plot([timeN(i), timeN(i+1)],[thetaDotDotSN(i), thetaDotDotSN(i+1)],'g');
    end
end
title('Angular Acceleration')
ylabel('Angular Acceleration (rad/s)')
xlabel('Time (s)')
xlim([0 timeN(end)])
ylim([-3000 3000]*pi/180)

% Plot Sled Movement
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

% Plot Experimental Phase Diagram
figure
for i = 1:size(moveIndex,2)-1
    hold on
    if moveIndex(i)==1
        plot([thetaN(i), thetaN(i+1)],[thetaDotSN(i), thetaDotSN(i+1)],'b');
    elseif moveIndex(i)==2        
        plot([thetaN(i), thetaN(i+1)],[thetaDotSN(i), thetaDotSN(i+1)],'r');
    else
        plot([thetaN(i), thetaN(i+1)],[thetaDotSN(i), thetaDotSN(i+1)],'g');
    end
end
title('Experimental Phase Diagram')
ylabel('Angular Velocity (rad/s)')
xlabel('Angular Displacement (rad)')
ylim([-1000 1000]*pi/180)


% Plot time lag embed plot
figure
for i = 1:size(moveIndex,2)-periodSteps-1
    hold on
    if moveIndex(i)==1
        plot([thetaN(i), thetaN(i+1)],[thetaN(i+periodSteps), thetaN(i+periodSteps+1)],'b');
    elseif moveIndex(i)==2
        plot([thetaN(i), thetaN(i+1)],[thetaN(i+periodSteps), thetaN(i+periodSteps+1)],'r');
    else
        plot([thetaN(i), thetaN(i+1)],[thetaN(i+periodSteps), thetaN(i+periodSteps+1)],'g');
    end
end
title('Time Lag Embed Diagram')
ylabel('Angular Displacement (rad)')
xlabel('Angular Displacement (rad)')

if storeVideo
    WriteVideo('Test2',double(storeVid))
end

if storeData
    save Damped2 timeN thetaN thetaDotSN thetaDotDotSN distN moveIndex
end
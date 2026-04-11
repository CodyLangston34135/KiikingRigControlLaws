function [storeThetaFit] = smoothData(theta,time,frameRate)

    [peaksP, locsP] = findpeaks(theta,'MinPeakProminence',50);
    [peaksN, locsN] = findpeaks(-theta,'MinPeakProminence',50);
    peaksN = -peaksN;

    plot(theta)
    hold on
    scatter(locsP,peaksP)
    scatter(locsN,peaksN)

    locs = [locsP, locsN];
    peaks = [peaksP, peaksN];

    [~, indLocs] = sort(locs,'ascend');

    for i = 1:length(indLocs)-1
        tempLoc(i,:) = [locs(indLocs(i)) locs(indLocs(i+1))];
        tempPeak(i,:) = [peaks(indLocs(i)) peaks(indLocs(i+1))];
    end

    locs = tempLoc;
    peaks = tempPeak;

    % figure
    storeThetaFit = zeros(1,1);
    for i = 1:size(locs,1)-1
        tempTheta = theta(locs(i,1):locs(i,2))';
        tempTime = time(locs(i,1):locs(i,2))'-time(locs(i,1));
        tempTimeN = tempTime(1):1/frameRate:tempTime(end);
        tempFreq = 1/(2*(time(locs(i,1))-time(locs(i,2))));
        
        % subplot(1,2,1)
        [tempTimeN,tempThetaFit] = testFit(tempTime,tempTheta,tempTimeN,false);

        % subplot(1,2,2)
        % plot(tempTimeN,tempThetaFit);
        % hold on
        % scatter(tempTime,tempTheta)
        
        storeThetaFit(1,end:end+length(tempThetaFit)-1) = tempThetaFit;
    end
end
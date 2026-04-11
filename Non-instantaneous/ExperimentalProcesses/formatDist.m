%% This is quite literally the worst way to do this

function [moveIndex, moveLength] = formatDist(distN,velN)

    tempIndex = abs(velN)>=0.18;
    for i = 1:size(distN,2)-3
        tempMat(i,1:4) = tempIndex(i:i+3);
    end
    tempIndex = ismember(tempMat,ones(1,4),'rows');
    
    moveIndex = zeros(1,size(distN,2));
    index = 1;
    moveLength = 0;
    for i = 1:length(tempIndex)
        if tempIndex(i) == 1 && velN(i) > 0
            moveIndex(1,i:i+3) = ones(1,4);
            moveLength(index) = moveLength(index)+1;
        elseif tempIndex(i) == 1 && velN(i) < 0
            moveIndex(1,i:i+3) = 2*ones(1,4);
            moveLength(index) = moveLength(index)+1;
        else
            index = index+1;
            moveLength(index) = 0;
        end
    end

    for i = 2:length(moveIndex)-1
        if moveIndex(i-1) == 0 && moveIndex(i) == 1 && moveIndex(i+1) == 2
            moveIndex(i) = 2;
        elseif moveIndex(i-1) == 2 && moveIndex(i) == 1 && moveIndex(i+1) == 0
            moveIndex(i) = 2;
        end
    end

    moveLength = nonzeros(moveLength)+3;
end
%% Plot Experimental
function plotExperimental(inputFile)

load(inputFile)

%% Plot Experimental Diagrams
% Plot experimental phase diagram
figure(1)
subplot(1,2,2)
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

% Plot Angular Displacement
figure(2)
subplot(2,2,2)
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
title('Experimental Angular Displacement')
ylabel('Angular Displacement (rad)')
xlabel('Time (s)')
xlim([0 timeN(end)])
ylim([-4 4])

% Plot Angular Velocity
subplot(2,2,4)
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
title('Experimental Angular Velocity')
ylabel('Angular Velocity (rad/s)')
xlabel('Time (s)')
xlim([0 timeN(end)])
ylim([-15 15])

% Plot time lag embed plot
periodSteps = 16;
figure(4)
subplot(1,2,2)
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
title('Experimental Time Lag Embed Diagram')
ylabel('Angular Displacement (rad)')
xlabel('Angular Displacement (rad)')
xlim([-4 4])
ylim([-4 4])

end
function [] = ...
    P300_EyeTracking_inMIClassifier(GazeX,GazeY, TrialTimers)

MonPos = get(groot,'MonitorPositions'); %Get monitor positions
Offset = 1000;  %Offset specified in EyetrackerLogger.cpp in ::OnGazeDataEvent
MonPos2 = MonPos(2,:); %Isolate second monitor
MonPos2([2 4]) = MonPos2([2 4])-282+Offset; %Yshift down 282 pixels

%% Define Gaze quality for each target
Targs = {'L','R','NG'};

%Measured target centers (in cm)
TargLocs = {[11.8 14.9],[36.1 14.9],[23.7 14.9]};
pxpercm = 35.5;
TargLocs = cellfun(@(x) x*pxpercm,TargLocs,'UniformOutput',false)

yoff = 1000-282;
xoff = 1366;


%%Plot

figure
scatter(GazeX,GazeY,[],[.8 .8 .8]); hold on;
line(MonPos2(1)*[1 1],[MonPos2(2) MonPos2(4)],'Color','k')
line((MonPos2(1)+MonPos2(3))*[1 1],[MonPos2(2) MonPos2(4)],'Color','k')
line([MonPos2(1) MonPos2(1)+MonPos2(3)],MonPos2(2)*[1 1],'Color','k')
line([MonPos2(1) MonPos2(1)+MonPos2(3)],MonPos2(4)*[1 1],'Color','k')
axis equal;
xlabel('X-postion'); ylabel('Y-position')
clrs = {[0 0 1],[0 1 0],[1 0 0]};
for j = 1:length(Targs)
    tmr = TrialTimers{j};
    targetbox = [xoff+TargLocs{j}(1)-50 yoff+TargLocs{j}(2)-50 ...
        xoff+TargLocs{j}(1)+50 yoff+TargLocs{j}(2)+50];
    line(targetbox(1)*[1 1],[targetbox(2) targetbox(4)],'Color',clrs{j});
    line(targetbox(3)*[1 1],[targetbox(2) targetbox(4)],'Color',clrs{j});
    line([targetbox(1) targetbox(3)],targetbox(2)*[1 1],'Color',clrs{j});
    line([targetbox(1) targetbox(3)],targetbox(4)*[1 1],'Color',clrs{j});
    text(targetbox(1),targetbox(2)+300,Targs{j},'Color',...
        clrs{j},'FontSize',20,'HorizontalAlignment','left',...
        'VerticalAlignment','bottom');
    for k = 1:size(tmr,1)    
        scatter(GazeX(tmr(k,1):tmr(k,2)),GazeY(tmr(k,1):tmr(k,2)),[],...
            (k/size(tmr,1))*clrs{j})
        gazeacc{j}(:,k) =  [mean(double(GazeX(tmr(k,1):tmr(k,2))))-...
            (xoff+TargLocs{j}(1));
            mean(double(GazeY(tmr(k,1):tmr(k,2))))-...
            (yoff+TargLocs{j}(2))];
        gazevar{j}(:,k) = [var(double(GazeX(tmr(k,1):tmr(k,2))));
            var(double(GazeY(tmr(k,1):tmr(k,2))))];
    end
end


figure
subplot(121); bar(cell2mat(cellfun(@(x) mean(x,2),gazeacc,'UniformOutput',false))');
title('Mean gaze offset');
ax = gca;
ax.XTickLabel = Targs;
subplot(122); bar(cell2mat(cellfun(@(x) mean(x,2),gazevar,'UniformOutput',false))');
title('Gaze variance');
ax = gca;
ax.XTickLabel = Targs;




end

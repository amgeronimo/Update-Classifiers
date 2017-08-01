function [gazeacc, gazevar, gazeinvalid, targetbox] = ...
    P300_EyeTracking_inP300Classifier_v2(GazeX,GazeY, Data, SC,SR, StimulusCode, ...
    SequencePhase, StimulusType, StateDuration, NumTrials, NumSequences, Speller, fs, c, ...
    Figures_On, Adopt_Monitor_Specs, TTSpell, AppPos)


%%
if Adopt_Monitor_Specs == 1 %Using Phillips Monitor
    MonPos =  [1           1        1920        1080;
        1921           1        1920        1080];
else
    MonPos = get(groot,'MonitorPositions'); %Get monitor positions
end

% AppPos = [c.WindowLeft.NumericValue c.WindowTop.NumericValue ...
%     c.WindowWidth.NumericValue c.WindowHeight.NumericValue];
if size(MonPos,1) == 2
    MonPos2 = MonPos(2,:); %Isolate second monitor
    ymondist = MonPos(2,4)-MonPos(1,4)+MonPos(2,2)-1; %find distance from top of monitor 1 to top of monitor 2
else
    MonPos2 = MonPos;
    ymondist = 0;
end

if (isfield(c,'LogEyetrackerTobiiX')||isfield(c,'LogEyeTrackerTobiiX'))
    COffset = 1000;  %Offset specified in EyetrackerLogger.cpp in ::OnGazeDataEvent
    AppPos(:,2) = AppPos(:,2)+COffset;
    MonPos2(2) = COffset-ymondist;
elseif (isfield(c,'LogEyeTribe') || isfield(c,'LogEyeTribe_v2'))
    mondiff = MonPos2(1:2)-[1 ymondist];
    AppPos(:,1:2) = AppPos(:,1:2)-mondiff;
    MonPos2(1) = 1; MonPos2(2) = 1;
elseif (isfield(c,'LogEyetrackerTobii3'))
       mondiff = MonPos2(1:2)-[1 ymondist];
    AppPos(:,1:2) = AppPos(:,1:2)-mondiff;
    MonPos2(1) = 1; MonPos2(2) = 1;
else
    disp('No Eyetrackers used.'); 
    gazeacc = []; gazevar = []; gazeinvalid = []; targetbox = [];
    return;
end


%Scale the input for Tobii X2-30 tracker
if isfield(c,'LogEyetrackerTobii3')
    GazeX = double(GazeX)/65535*MonPos2(3) + MonPos2(1);
    GazeY = double(GazeY)/65535*MonPos2(4) + MonPos2(2);
end

% figure
% subplot(221)
% scatter(1:length(GazeX),GazeX,[],parula(length(GazeX)))
% line([1 length(GazeX)],MonPos2(1)*[1 1],'Color','b')
% text(length(GazeX)/2,MonPos2(1)*.9,'Left','HorizontalAlignment','center')
% line([1 length(GazeX)],[MonPos2(1)+MonPos2(3)]*[1 1],'Color','r')
% text(length(GazeX)/2,MonPos2(1)+MonPos2(3)*1.1,'Right','HorizontalAlignment','center')
% ylabel('X-postion');
% subplot(223)
% scatter(1:length(GazeY),GazeY,[],parula(length(GazeY)))
% line([1 length(GazeY)],MonPos2(2)*[1 1],'Color','g')
% text(length(GazeY)/2,MonPos2(2)*1.1,'Top','HorizontalAlignment','center')
% line([1 length(GazeY)],(MonPos2(2)+MonPos2(4))*[1 1],'Color','k')
% text(length(GazeY)/2,(MonPos2(2)+MonPos2(4))*1.1,'Bottom','HorizontalAlignment','center')
% ylabel('Y-position');
% subplot(2,2,2)
% scatter(GazeX,GazeY,[],parula(length(GazeX)))
% line(MonPos2(1)*[1 1],[MonPos2(2) MonPos2(2)+MonPos2(4)],'Color','b')
% line((MonPos2(1)+MonPos2(3))*[1 1],[MonPos2(2) MonPos2(2)+MonPos2(4)],'Color','r')
% line([MonPos2(1) MonPos2(1)+MonPos2(3)],MonPos2(2)*[1 1],'Color','g')
% line([MonPos2(1) MonPos2(1)+MonPos2(3)],(MonPos2(2)+MonPos2(4))*[1 1],'Color','k')
% axis equal;
% xlabel('X-postion'); ylabel('Y-position')




%% Define Gaze quality for each letter

%Define Beginning and End of trials
tmp = find(SequencePhase == 2);
tmp2 = diff(tmp);
tmp3 = find(tmp2>1);
trialboundaries = [[tmp(1); tmp(tmp3+1)] [tmp(tmp3); tmp(end)]];

%Define targets for each trial
targetstim = double(StimulusCode(logical(StimulusType)));
targetstim = targetstim(1:StateDuration:length(targetstim));
 k = 1;
 ki = 1;
 for r = 1:length(NumTrials)
    for i = 1:NumTrials{r}
        ts(ki,:) = unique(targetstim(k:k+NumSequences{r}*2-1))';
        k=k+NumSequences{r}*2;
        ki=ki+1;
    end
end
targetstim = ts;

TargH = c.TargetHeight.NumericValue;
TargW = c.TargetWidth.NumericValue;
BarH = c.StatusBarSize.NumericValue;
NumR = c.NumMatrixRows.NumericValue;
NumC = c.NumMatrixColumns.NumericValue;

%Find locations of these targets
if Speller == 1 %Checkerboard speller
    %Map the target codes to the stimulus
    clear lettertarg
    for ii = 1:length(SC)
    map = [cellfun(@(x) str2num(x),SR{ii}) cellfun(@(x) str2num(x),SC{ii})];
    lettertarg(ii) = find(sum(targetstim(ii,2)==map | targetstim(ii,1)==map,2)>1);
    end
    
    %Didnt want to figure out how to automatically generate the locations,
    %so did it manually
%     letterlocs = [1 2 3 4 5 6 1 2 3 4 5 6 1 2 5 6 1 2 5 6 1 2 3 4 5 6 1 2 3 4 5 6;
%         1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6];
    letterlocs = [1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8;
                  1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5];
    targlocs = letterlocs(:,lettertarg)';

    %This works for the case in which there is a 6x6 grid, with 32 spots
    %occupied by icons, with target width =6.25% and vertical
    %target height = 10%.  This is 1.6 ratio, same as the monitor.
    %For the monitor, there are 35.5 pixels per cm.
%     pxpercm = 35.5;
%     targleft = 14.8*pxpercm;
%     targtop = 7.05*pxpercm;
%     cmpertarg = 2.95; %The width(and height) of targets.

%   targleft = 15.9*pxpercm;
%      targtop = 9.2*pxpercm;
%     cmpertarg = 2.58; %The width(and height) of targets. Found this by
%     %measuring the height and width of the block (15.5cm) and dividing by 6

    temp_dist_w = ((100-(TargW*NumC))/2)/100*AppPos(1,3);
    temp_dist_h = (((100-BarH)-(TargH*NumR))/2)/100*AppPos(1,4)+BarH/100*AppPos(1,4);

    
elseif Speller == 2 %Single Flash (RSVP) speller.
    lettertarg = min(targetstim,[],2);
    targlocs = ones(length(lettertarg),2);
    
    %These are the dimensions in the case of the RSVP speller, using
    %6.25% width and 10% height
    %     pxpercm = 35.5;
    
    temp_dist_w = ((100-TargW)/2)/100*AppPos(1,3);
    
    temp_dist_h = (((100-BarH)-TargH)/2)/100*AppPos(4)+BarH/100*AppPos(1,4);
    
elseif Speller == 3 %Grid Speller (default P3 task)
    temp = sort(targetstim,2);
    temp(:,2) = temp(:,2)-6;
    lettertarg = 6*(temp(:,1)-1)+temp(:,2);
    targlocs = fliplr(temp);
    
    %       pxpercm = 35.5;
    %     targleft = 14.8*pxpercm;
    %     targtop = 7.05*pxpercm;
    %     cmpertarg = 2.95; %The width(and height) of targets.
    temp_dist = ((100-TargSize)/2)/100*AppPos(1,4);

end

targleft = temp_dist_w+AppPos(:,1);
targtop = temp_dist_h+AppPos(:,2);
pxtarg_h = TargH/100*AppPos(:,4);
pxtarg_w = TargW/100*AppPos(:,3);

%%Plot
TrlInd = regexp(TTSpell,'[A-Z,1-8]');
TargSymb = cell(length(TrlInd),1);
for kk = 1:length(TrlInd)
    if kk == length(TrlInd)
        TargSymb{kk} = TTSpell(TrlInd(kk):end);
    else
        TargSymb{kk} = TTSpell(TrlInd(kk):TrlInd(kk+1)-1);
    end
end


%Invert GazeY data and the Y Monitor positions (also invert target Y
%positions in the loop below)
GazeY = -double(GazeY);
MonPos2([2 4]) = -MonPos2([2 4]);
AppPos(:,[2 4]) = -AppPos(:,[2 4]);

if Figures_On == 1
figure('Name','2D position of targets and gaze')
scatter(GazeX,GazeY,[],[.8 .8 .8],'.'); hold on;
line(MonPos2(1)*[1 1],[MonPos2(2) MonPos2(2)+MonPos2(4)],'Color','k')
line((MonPos2(1)+MonPos2(3))*[1 1],[MonPos2(2) MonPos2(2)+MonPos2(4)],'Color','k')
line([MonPos2(1) MonPos2(1)+MonPos2(3)],MonPos2(2)*[1 1],'Color','k')
line([MonPos2(1) MonPos2(1)+MonPos2(3)],(MonPos2(2)+MonPos2(4))*[1 1],'Color','k')
line(AppPos(1,1)*[1 1],[AppPos(1,2) AppPos(1,2)+AppPos(1,4)],'Color',[.5 .5 .5])
line((AppPos(1,1)+AppPos(1,3))*[1 1],[AppPos(1,2) AppPos(1,2)+AppPos(1,4)],'Color',[.5 .5 .5])
line([AppPos(1,1) AppPos(1,1)+AppPos(1,3)],AppPos(1,2)*[1 1],'Color',[.5 .5 .5])
line([AppPos(1,1) AppPos(1,1)+AppPos(1,3)],(AppPos(1,2)+AppPos(1,4))*[1 1],'Color',[.5 .5 .5])
line([AppPos(1,1) AppPos(1,1)+AppPos(1,3)],(AppPos(1,2)+BarH/100*AppPos(1,4))*[1 1],'Color',[.5 .5 .5])
axis equal;
xlabel('X-postion'); ylabel('Y-position')
clrs = parula(length(lettertarg));

for j = 1:length(lettertarg)

    templeft = targleft(j) + (targlocs(j,1)-1)*pxtarg_w(j);
    temptop = targtop(j) + (targlocs(j,2)-1)*pxtarg_h(j);
    targetbox = [templeft -temptop templeft+pxtarg_w(j)...
         -temptop-pxtarg_h(j)];
     line(targetbox(1)*[1 1],[targetbox(2) targetbox(4)],'Color',clrs(j,:));
line(targetbox(3)*[1 1],[targetbox(2) targetbox(4)],'Color',clrs(j,:));
line([targetbox(1) targetbox(3)],targetbox(2)*[1 1],'Color',clrs(j,:));
line([targetbox(1) targetbox(3)],targetbox(4)*[1 1],'Color',clrs(j,:));
scatter(GazeX(trialboundaries(j,1):trialboundaries(j,2)),...
    GazeY(trialboundaries(j,1):trialboundaries(j,2)),'.',...
    'MarkerEdgeColor',clrs(j,:));
text(targetbox(1),targetbox(2),TargSymb{j},'Color',...
    clrs(j,:),'FontSize',18,'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
end
end


%% Plot eye movements vs EOG channels
gX = double(GazeX);

if Figures_On == 1
figure('Name','Eye position with EOG channels'); hold on;
plot(Data');
% plot(diff(Data));
plot(100*(gX-mean(gX))/range(gX),'k','LineWidth',2);
end

%% Plot the periods of eye invalidity

DgX = gX-circshift(gX,1);
%find stretches of diff(gX) that are at least .5 seconds long
DgX2 = DgX'~=0;
DgX2 = [1 DgX2 1];
frontSTR = strfind(DgX2,[1 0]);
backSTR = strfind(DgX2,[0 1]);

repeat0length = backSTR-frontSTR;
abovelimit = find(repeat0length>(fs/2));

invaliddata = zeros(length(DgX2),1);
for i = abovelimit
    %the -1 is to remove the 1 added earlier in DgX2
    invaliddata(frontSTR(i):frontSTR(i)+repeat0length(i)-1)=1;
end

if Figures_On == 1
figure('Name','Gaze X position, its derivative, and time of eye gaze invalidity'); hold on;
plot(gX);
plot(DgX);
plot(invaliddata*5000);
end

%Calculate the percentage of the trial for which the gaze was valid
gazeinvalid = zeros(length(lettertarg),1);
for j = 1:length(lettertarg)
gazeinvalid(j) =  sum(invaliddata(trialboundaries(j,1):trialboundaries(j,2)))...
    /diff(trialboundaries(j,:))*100;
end

%% Gaze Accuracy and Variance
gazeacc = zeros(length(lettertarg),2);
gazevar = zeros(length(lettertarg),2);
targetbox = zeros(length(lettertarg),4);
for j = 1:length(lettertarg)
    templeft = targleft(j) + (targlocs(j,1)-1)*pxtarg_w(j);
    temptop = targtop(j) + (targlocs(j,2)-1)*pxtarg_h(j);
    targetbox(j,:) = [templeft -temptop templeft+pxtarg_w(j)...
         -temptop-pxtarg_h(j)];
    gxdata = double(GazeX(trialboundaries(j,1):trialboundaries(j,2)));
    gxdata(logical(invaliddata(trialboundaries(j,1):trialboundaries(j,2)))) = NaN;
    gydata = double(GazeY(trialboundaries(j,1):trialboundaries(j,2)));
    gydata(logical(invaliddata(trialboundaries(j,1):trialboundaries(j,2)))) = NaN;
    gazeacc(j,:) =  [nanmean(gxdata)-(templeft+.5*pxtarg_w(j));
        nanmean(gydata)-(-temptop-.5*pxtarg_h(j))];
    gazevar(j,:) = [nanvar(gxdata); nanvar(gydata)];
end
%%
if Figures_On == 1
    figure('Name','Gaze Statistics','Position',[100 100 900 400])
    subplot(1,3,1); bar(gazeacc); title('Mean gaze offset');
    ax = gca; yll = ax.YLim; ax.XTickLabel = [];
    for i = 1:length(lettertarg)
        text(i,yll(1)-diff(yll)/30,TargSymb{i},...
            'HorizontalAlignment','center');
    end
    subplot(1,3,2); bar(gazevar); title('Gaze variance');
    ax = gca; yll = ax.YLim; ax.XTickLabel = [];
    for i = 1:length(lettertarg)
        text(i,yll(1)-diff(yll)/30,TargSymb{i},...
            'HorizontalAlignment','center');
    end
    subplot(1,3,3); bar(gazeinvalid); title('Percent invalid gaze');
    ax = gca;
    ax = gca; yll = ax.YLim; ax.XTickLabel = [];
    for i = 1:length(lettertarg)
        text(i,yll(1)-diff(yll)/30,TargSymb{i},...
            'HorizontalAlignment','center');
    end
end


end

%Files updated June/July 2015 for use with second round of ALS patients
%In addition to some minor changes, added the ability to use a regularized
%(LASSO) regression classifier option (Cfier == 4).
function [CM, acc_fromCfier] = BuildMIClassifier_teleBCI(Name, Session, Run, BuildC,Cfier,...
    Batch,RerefVal,Art,DoCSP,WhichDiff,CSession,CRun)
close all

%% Parameters

%%MAJOR changes to file since beginning of recordings in ALS clinic.


%11/15/13 - changed the window on Pfeat to use start+.5:
%end-1.5 instead of start+1:end-1.  This matches up with the windows used
%during feedback, because tT starts at .5!  All Batch_run results prior to
%this date will be different for this reason.
cdp = cd;
[~, sl] = regexp(cdp,'\Users\');
se = regexp(cdp,'\');
tmp = find(sl==se)
pcusr = cdp(sl+1:se(tmp+1)-1);
if strcmp(pcusr,'ageronimo')
    bcifold = '\Documents\BCI2000_305_VS2012';
else
    bcifold = '\My Documents\Year4\BCI2000_305_2015';
end
NameBase = ['C:\Documents and Settings\' pcusr bcifold];

Figures_On = 0;
AvgF = 0;

if Batch == 0
    WhichDiff = 1; %1 = L-R, 2 = L-C, 3 = R-C
    RerefVal = 2; %(0 - no reref, 1 - CAR, 2 - LAP)
    Art=2; %Artifact Data (0 - nothing, 1 - remove regions of artifact, 2 - regression)
    DoCSP = 0;
    Figures_On = 1;
end

% EOGloc  = [];
% EEGloc = 1:2;
% EOGloc  = 15:16;
% EEGloc = 1:14;
% EOGloc = []
% EEGloc = 1:5;
EOGloc = 20:22;
EEGloc = 1:19;
% EOGloc = [];


if BuildC == 0
    if Batch == 0
        CSession = input('What Session? (input as string)   ');
        CRun = input('What Run? (input as string)   ');
    elseif Batch == 1
        CRun = cell2mat(CRun);
    end
    %     NameBase = ['C:\Documents and Settings\amg5106\My Documents\'...
    %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' Name CSession '\'...
    %         Name 'S' CSession 'R' CRun '_'];
    NameBaseUserData = [NameBase '\data\' ...
        Name CSession '\' Name 'S' CSession 'R' CRun '_'];
end

%% Load Data
tic
Data = [];
Cursor = [];
ResCode = [];
TC = [];
Feedback = [];
NewRunName = [];
GazeX = [];
GazeY = [];
NT = 0;
for rr = 1:length(Run)
    %     [a b c d] = load_bcidat(['C:\Documents and Settings\amg5106\My Documents\'...
    %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' ...
    %         Name Session '\' Name 'S' Session 'R' Run{rr} '.dat'],'-calibrated');
    try
    [a b c d] = load_bcidat([NameBase '\data\' ...
        Name Session '\' Name 'S' Session 'R' Run{rr} '.dat'],'-calibrated');
    disp(['Time to load data ' num2str(toc)]);
    catch
    [a b c d] = load_bcidat(['P:\ALS Proj Data\' ...
        Name Session '\' Name 'S' Session 'R' Run{rr} '.dat'],'-calibrated');
    disp(['Time to load data ' num2str(toc)]);    
    end
    
    % Shift Data? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         kk = 1;
%         for jj = [1 21]
%             [MM, ~] = peakdet(a(:,jj),100);
%             MaxTab{kk} = MM(:,1);
%             kk = kk+1;
%         end
%     
%         %Positive differences are when amp2 leads amp1.
%         CumDiff = []; CumLoBc = [];
%         for jj = 1:length(MaxTab{1})
%             CumDiff = [CumDiff; MaxTab{1}(jj)-MaxTab{2}];
%             CumLoc = [CumLoc; repmat(MaxTab{1}(jj),length(MaxTab{2}),1)];
%         end
%     
%         keepem = CumDiff<10&CumDiff>-10;
%         CumDiff = CumDiff(keepem);
%         CumLoc = CumLoc(keepem);
%     
%     
%         if Figures_On == 1
%             figure
%             tx(1) = subplot(221); plot(a(:,1)); hold on; plot(a(:,21),'r');
%             subplot(2,2,[2 4]); hist(CumDiff,[-10:10]); xlim([-10 10])
%         end
%     
%         Mshift = mode(CumDiff);
%         if Mshift > -2 && Mshift < 2
%         a(:,17:end) = circshift(a(:,17:end),[Mshift 0]);
%         disp(['THE DATA HAS BEEN SHIFTED ' num2str(Mshift)])
%         end
%     
%         if Figures_On == 1
%             tx(2) = subplot(223); plot(a(:,1)); hold on; plot(a(:,21),'r');
%             title(['Shift ' num2str(Mshift)]);
%             linkaxes(tx,'x')
%         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fs = c.SamplingRate.NumericValue;
    sbs = c.SampleBlockSize.NumericValue;
    if BuildC == 0
        %For replication of online values, dont cut off any of the data
        cutoff = 1;
        Data = [Data;  a(cutoff:end,:)];
        Cursor = [Cursor; double(b.CursorPosX(cutoff:end))];
        ResCode = [ResCode; double(b.ResultCode(cutoff:end))];
        tempTC = double(b.TargetCode(cutoff:end)); tempTC(1:c.PreRunDuration.NumericValue*fs+9)=0;
        TC = [TC; tempTC];
        Feedback = [Feedback; double(b.Feedback(cutoff:end))];
        if isfield(b,'EyetrackerLeftEyeGazeX')
            GazeX = [GazeX; b.EyetrackerLeftEyeGazeX(cutoff:end)];
            GazeY = [GazeY; b.EyetrackerLeftEyeGazeY(cutoff:end)];
        end
        
    elseif BuildC == 1
        %Trim Data and StimulusCodes if building classifier
        cutoff = c.PreRunDuration.NumericValue*fs+1;
        Data = [Data;  a(cutoff:end,:)];
        Cursor = [Cursor; double(b.CursorPosX(cutoff:end))];
        ResCode = [ResCode; double(b.ResultCode(cutoff:end))];
        tempTC = double(b.TargetCode(cutoff:end)); tempTC(1:sbs)=0;  %Overwrite the first points just in case there is overflow from the previous run.
        TC = [TC; tempTC];
        Feedback = [Feedback; double(b.Feedback(cutoff:end))];
        if isfield(b,'EyetrackerLeftEyeGazeX')
            GazeX = [GazeX; b.EyetrackerLeftEyeGazeX(cutoff:end)];
            GazeY = [GazeY; b.EyetrackerLeftEyeGazeY(cutoff:end)];
        end
    end
    
    NewRunName = strcat(NewRunName,Run{rr});
    
    %Load Classifier
    ClassifierUsed{rr} = c.Classifier.NumericValue;
    ClassifierUsed{rr}(:,4) = cellfun(@str2num,c.Classifier.Value(:,4));
    
    %Number of Channels
    NC{rr} = size(Data,2);
    %Number of Trials
    NT = NT + c.NumberOfTrials.NumericValue;
    
end
Data = Data';
%% Check Variables

%Are classifiers/channels consistant across runs?
if sum(diff(cell2mat(NC))) ~= 0
    disp('Different number of channels in each run!');
end
NC = NC{1};

if NC ~= max([EEGloc EOGloc])
    disp('Channels do not match specified EEG and EOG locations');
end


try %Are classifiers the same size
    c2mC = cell2mat(ClassifierUsed);
    if sum(sum(diff(reshape(c2mC,size(c2mC,1),size(c2mC,2)/rr,rr),[],3))) ~= 0
        disp('Classifiers have different values in each run!');
    end
catch err
    disp('Classifiers are different sizes in each run!');
end
%If the same size, do they contain the same values?

ClassifierUsed = ClassifierUsed{1};


%% Find superlarge (>300 mV) artifacts

%If not using the contents of this cell, set GoodData to ones;
GoodData = ones(size(Data,2),1);


% BadData = [];
% for i = 1:size(Data,1)
%     %Find large artifacts
%     BadData = [BadData find(Data(i,:)>300|Data(i,:)<-300)];
% end
%     %Define a range 1 second before and 3 seconds after all of the found
%     %artifacts
%     BadData = repmat(BadData,4*fs,1)+repmat((-1*fs+1:3*fs)',1,length(BadData));
%     BadData = unique(BadData);
%     GoodData = ~ismember(1:length(Data),BadData);


%% Plot raw data and trim bad data

% if Figures_On == 1
%     figure
%     for i = 1:size(Data,1)
%         if size(Data,1) == 22
%             eleclocs = [7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 29 2 3 4];
%             tr(i) = subplot(6,5,eleclocs(i)); plot(Data(i,:));
%             hold on;
%             plot(100*(1-GoodData),'r');
%         else
%             tr(i) = subplot(6,5,i); plot(Data(i,:));
%             hold on;
%             plot(100*(1-GoodData),'r');
%         end
%     end
%     linkaxes(tr);
% end
% 
% GoodData = logical(GoodData);
% 
% Data = Data(:,GoodData);
% Cursor = Cursor(GoodData);
% ResCode = ResCode(GoodData);
% TC = TC(GoodData);
% Feedback = Feedback(GoodData);


%% Triggers

%Find locations of left and right cues
LC = TC==2;
RC = TC==1;
CC = TC==3;

if Figures_On==1
    figure
    mm(1)=subplot(311);plot(LC);
    mm(2)=subplot(312);plot(RC,'r');
    mm(3)=subplot(313);plot(CC,'g');
    linkaxes(mm)
end

% if BuildC == 0
%     LC(1:c.PreRunDuration.NumericValue*fs)=0;
%     RC(1:c.PreRunDuration.NumericValue*fs)=0;
%     CC(1:c.PreRunDuration.NumericValue*fs)=0;
% elseif BuildC == 1
%     LC(1:9)=0;
%     RC(1:9)=0;
%     CC(1:9)=0;
% end




LCt = LC-circshift(LC,[1 0]);
RCt = RC-circshift(RC,[1 0]);
CCt = CC-circshift(CC,[1 0]);


%Set the first value to zero, since one of them will have a -1 in the first
%slot
if RCt(1)==-1
    RCt = circshift(RCt,[-1 0]);
end
if LCt(1)==-1
    LCt = circshift(LCt,[-1 0]);
end
if CCt(1)==-1
    CCt = circshift(CCt,[-1 0]);
end



if Figures_On == 1
    figure
    mm(1)=subplot(311);plot(LCt);
    mm(2)=subplot(312);plot(RCt,'r');
    mm(3)=subplot(313);plot(CCt,'g');
    linkaxes(mm)
end


L_trialstart = find(LCt==1);
L_trialend = find(LCt==-1);
R_trialstart = find(RCt==1);
R_trialend = find(RCt==-1);
C_trialstart = find(CCt==1);
C_trialend = find(CCt==-1);
%Removing unfinished trials if necessary

if BuildC == 1 %added 11/1/13 -- only need to make sure there is enough data when building classifier.
incL = find(L_trialend-L_trialstart<768); %make sure there is more than three seconds of trial
L_trialstart = L_trialstart(~ismember(1:length(L_trialstart),incL));
L_trialend = L_trialend(~ismember(1:length(L_trialend),incL));

incR = find(R_trialend-R_trialstart<768);
R_trialstart = R_trialstart(~ismember(1:length(R_trialstart),incR));
R_trialend = R_trialend(~ismember(1:length(R_trialend),incR));

incC = find(C_trialend-C_trialstart<768);
C_trialstart = C_trialstart(~ismember(1:length(C_trialstart),incC));
C_trialend = C_trialend(~ismember(1:length(C_trialend),incC));
end

% if length(L_trialstart)==1+length(L_trialend)
%     L_trialstart = L_trialstart(1:end-1);
% elseif length(L_trialstart)>1+length(L_trialend)
%     disp('problem with trial markers'); return;
% else
% end
% if length(R_trialstart)==1+length(R_trialend)
%     R_trialstart = R_trialstart(1:end-1);
% elseif length(R_trialstart)>1+length(R_trialend)
%     disp('problem with trial markers'); return;
% else
% end
% if length(C_trialstart)==1+length(C_trialend)
%     C_trialstart = C_trialstart(1:end-1);
% elseif length(C_trialstart)>1+length(C_trialend)
%     disp('problem with trial markers'); return;
% else
% end

%Equating left, right, and no cue trials if necessary -- REMOVING THIS
%BECAUSE IT LIMITS THE TRIALS WE CAN USE WHEN WE USE THIS FUNCTION FOR
%MULTIPLE RUNS.
NTl = length(L_trialstart);
NTr = length(R_trialstart);
NTc = length(C_trialstart);
% if NTc ~= 0
%     [minz, iminz] = min([NTl NTr NTc]);
%     L_trialstart = L_trialstart(1:minz); L_trialend = L_trialend(1:minz); NTl = minz;
%     R_trialstart = R_trialstart(1:minz); R_trialend = R_trialend(1:minz); NTr = minz;
%     C_trialstart = C_trialstart(1:minz); C_trialend = C_trialend(1:minz); NTc = minz;
% else
%     [minz, iminz] = min([NTl NTr]);
%     L_trialstart = L_trialstart(1:minz); L_trialend = L_trialend(1:minz); NTl = minz;
%     R_trialstart = R_trialstart(1:minz); R_trialend = R_trialend(1:minz); NTr = minz;
% end


%End function if the number of trial types does not match the WhichDiff
%parameters
switch WhichDiff
    case 1
        if NTl==0 || NTr == 0
            CM = NaN; acc = NaN;
            return
        end
    case 2
        if NTl==0 || NTc == 0
            CM = NaN; acc = NaN;
            return
        end
    case 3
        if NTc==0 || NTr == 0
            CM = NaN; acc = NaN;
            return
        end
end




%% Calculate ERPs
% %If you get an error here, it is likely because the trial didnt finish and
% %this portion of code is trying to access data 2 seconds after the start of
% %a trial that didnt finish.
% if BuildC == 1
%     for li = 1:length(L_trialstart)
%         LData(li,:,:) = Data(:,L_trialstart(li):L_trialstart(li)+2*fs-1);
%     end
%     for ri = 1:length(R_trialstart)
%         RData(ri,:,:) = Data(:,R_trialstart(ri):R_trialstart(ri)+2*fs-1);
%     end
%     for ci = 1:length(C_trialstart)
%         CData(ci,:,:) = Data(:,C_trialstart(ci):C_trialstart(ci)+2*fs-1);
%     end
%     figure
%     pi = 1;
%     if size(Data,1)==16
%         VEPchans = [1 2 5 10 13 14 15 16];
%         VEPlabels = {'F3','F4','Cz','Pz','O1','O2','lEOG','cEOG'};
%     elseif size(Data,1)==22
%         VEPchans = [1 2 10 15 18 19 20 21 22];
%         VEPlabels = {'FP1','FP2','Cz','Pz','O1','O2','lEOG','cEOG','rEOG'};
%     end
%
%     %Target appears at 0 seconds, then cursor at 1 second, feedback runs from
%     %1-4 seconds, and then both target and cursor leave at second 5.
%     for ci = VEPchans
%         subplot(3,3,pi); plot((1:2*fs)/fs,squeeze(mean(LData(:,ci,:),1))); hold on;
%         plot((1:2*fs)/fs,squeeze(mean(RData(:,ci,:),1)),'r');
%         plot((1:2*fs)/fs,squeeze(mean(CData(:,ci,:),1)),'g');
%         title(VEPlabels{pi});
%         pi = pi+1;
%     end
%
% end


%% Artifact Reduction
if Art == 1 % Locate times with artifact and automatically remove
    
    
    
    %The oreder of indices may need to switch on all the references to Data
    try
        ArtInt = find(Data(:,EOGloc(1))>70&circshift(Data(:,EOGloc(1)),[1 0])<70);
    catch
        disp('Cannot remove data because no EOG channels were recorded'); CM =[];
        return
    end
    Arange = fix(-.1*fs:1*fs);
    ArtInt = repmat(ArtInt,1,length(Arange))+repmat(Arange,length(ArtInt),1);
    ArtInt = ArtInt';
    ArtInt = ArtInt(:);
    ArtInt = ArtInt(ArtInt>0 & ArtInt<size(Data,1));
    ArtLoc = zeros(size(Data,1),1);
    ArtLoc(ArtInt)=1;
    
    if Figures_On == 1
        figure
        plot(Data(:,EOGloc(1))); hold on; plot(50*ArtLoc,'r'); plot(51*TC,'g');
    end
    
    %This is the data with the Artifacts removed
    Data = Data(~ArtLoc,:);
    %Do the same for other states as well (Taking care of triggers below)
    Cursor = Cursor(~ArtLoc);
    ResCode = ResCode(~ArtLoc);
    Feedback = Feedback(~ArtLoc);
    TC = TC(~ArtLoc);
    
    if Figures_On == 1
        figure
        plot(Data(:,EOGloc(1))); hold on; plot(50*TC(~ArtLoc),'g');
    end
    
    %Now shift the triggers back
    for sb = 1:length(L_trialstart)
        moveback = sum(ArtLoc(1:L_trialstart(sb)));
        L_trialstart(sb) = L_trialstart(sb)-moveback+1;
        moveback = sum(ArtLoc(1:L_trialend(sb)));
        L_trialend(sb) = L_trialend(sb)-moveback;
        moveback = sum(ArtLoc(1:R_trialstart(sb)));
        R_trialstart(sb) = R_trialstart(sb)-moveback+1;
        moveback = sum(ArtLoc(1:R_trialend(sb)));
        R_trialend(sb) = R_trialend(sb)-moveback;
        scatter(L_trialstart(sb),50,'bo');
        scatter(L_trialend(sb),50,'bx');
        scatter(R_trialstart(sb),100,'ro');
        scatter(R_trialend(sb),100,'rx');
    end
    
    ArtWeights = eye(NC);
    ArtWeights2 = eye(NC,length(EEGloc));
    NC = length(EEGloc);
    
    
elseif Art == 2 %Artifact regression in the line of Schlogl et al.
    disp('Artifact Reduction')
    if BuildC == 1
        %Build Artifact rejector (output chans x input chans)
        %Find the eog & eeg channels
        
        %     EOGc = regexp(c.ChannelNames.Value,'EOG');
        %     Cloc = cellfun(@isempty,EOGc);
        %     EEGloc = find(Cloc~=0);
        %     EOGloc = find(Cloc==0);
        %
        %         %use the regression algorithm
        %     EOGD = Data(1:5*fs,:);
        
        %         figure
        %         plot(Data(:,EOGloc(1)));
        %         EOGDloc = input('Data points to use for EOG regression: ');
        
        try
            ArtInt = find(Data(EOGloc(1),:)-Data(EOGloc(2),:)>75 | ...
                Data(EOGloc(1),:)-Data(EOGloc(2),:)<-75);
        catch
            disp('Cannot remove data because no EOG channels were recorded');
            return
        end
        
        Arange = fix(-.1*fs:.5*fs);
        
        ArtInt = repmat(ArtInt',1,length(Arange))+repmat(Arange,length(ArtInt),1);
        ArtInt = ArtInt';
        ArtInt = unique(ArtInt(:));
        ArtInt = ArtInt(ArtInt>0 & ArtInt<size(Data,2));
        ArtLoc = logical(zeros(size(Data,2),1));
        ArtLoc(ArtInt) = true;
        
        if Figures_On == 1
            figure
            plot(Data(EOGloc(1),:)-Data(EOGloc(2),:))
            hold on
            plot(100*ArtLoc,'r')
        end
        
        EOGD = Data(:,ArtLoc);
        
        %the function covm adds an additional column of ones in front of the data
        %and is necessary for regress_eog.m
        if length(EOGloc)==3 && ~isempty(EOGD)
            [R] = regress_eog(covm(EOGD','E'),EEGloc, ...
                sparse([EOGloc(2),EOGloc(1),EOGloc(2),EOGloc(3)],[1,1,2,2],[1,-1,1,-1]));
        elseif length(EOGloc)==2 && ~isempty(EOGD)
            [R] = regress_eog(covm(EOGD','E'),EEGloc, ...
                sparse([EOGloc(1),EOGloc(2)],[1,1],[1,-1]));
        else
            R.r0 = eye(size(Data,1));
        end
        %Create full matrix for online artifact reduction
        %I believe this is the way they say to do it (pad Data with a channel
        %of ones -- this introduces a bias to the output channel) (see DD2 below).
        %However, this padding is not something I want to do online, and since
        %it is only a bias, we can remove the first column of ArtWeights.
        ArtWeights = full(R.r0)';
        ArtWeights2 = ArtWeights(EEGloc,:);
        NC = length(EEGloc);
        %ArtWeights = full(R.r1);
        %DD2 = [ones(size(Data,1),1),Data] * ArtWeights';
        
    elseif BuildC == 0
        
        
    end
    
    %Using the correction coefficients, transform entire training run to reduce
    %artifact
    
else
    ArtWeights = eye(size(Data,1));
    ArtWeights2 = eye(length(EEGloc),size(Data,1));
    NC = length(EEGloc);
end

if Figures_On==1
    if BuildC == 1
        Dw = ArtWeights2*Data;
        %
        %
        %Plot the F and O raw and artifacted, as well as EOG channels
        figure
        tt(1)=subplot(311); plot(Data(1,:)); hold on; plot(Dw(1,:),'r');
        tt(2)=subplot(312); plot(Data(2,:)); hold on;
        plot(Dw(2,:),'r');
        tt(3)=subplot(313); plot(Data(EOGloc,:)'); hold on
        linkaxes(tt,'x');
        
        clear Dw;
    end
end


%% Rereference Data
switch RerefVal
    case 0
        RefFilt = eye(length(EEGloc));
        %         RefFilt = eye(NC);
        %         DataF = Data*RefFilt;
    case 1
        RefFilt = (-1/length(EEGloc))*(ones(length(EEGloc))-eye(length(EEGloc)));
        RefFilt = RefFilt + eye(length(EEGloc));
        %         RefFilt = blkdiag(RefFilt, eye(length(EOGloc)));
        %         DataF = Data*RefFilt;
        
    case 2
        %For this laplacian configuration, Channels are F3, F4, T7, C3, Cz,
        %C4, T8, P7, P3 Pz, P4, P8, O1 ,O2.  Outputs are C3,C4,P3,P4.
        if size(Data,1)==16
            RefFilt = [-1/4 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0 0 0;
                0 -1/4 0 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0;
                0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 0 -1/4 0;
                0 0 0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 -1/4];
        elseif size(Data,1)==22
            %This is for the 22 channel setup
            RefFilt = [-1/4 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0 0 0 0 0 0 0 0;
                0 -1/4 0 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0 0 0 0 0 0;
                0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0 0 0;
                0 0 0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 0 0 -1/4 0 0 0;
                0 0 0 0 0 0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 0 -1/4 0;
                0 0 0 0 0 0 0 0 0 0 -1/4 0 0 0 -1/4 1 -1/4 0 -1/4];
        end
        NC_old = NC;
        NC = size(RefFilt,1);
end


% figure
% for ch = 1:NC
%     mmd(ch) = subplot(4,4,ch); plot(Data(:,ch));
%     hold on; plot(DataF(:,ch),'r');
% end
% linkaxes(mmd,'x');

%% Finding Reactive Band

% %I DONT BELIEVE ANYTHING USEFUL HAPPENS HERE - i have no way (at this
% %moment) of implementing a filter in BCI2000) -- Since the features to the
% %classifier are the instantaneous power in each frequency band, there is no
% %benefit of prefiltering in a reactive band.
%
Frange = [6 30];
%
% if BuildC == 1
%     %Data window is 128 points long, multiplied by (half of) a hamming window
%     wndw = c.FFTWindowLength.NumericValue*fs;
%     hm = hamming(wndw*2+1);
%     hm = hm(wndw+1:end-1);
%     Data_R = RefFilt*ArtWeights2*Data;
%     PL = zeros(NC,wndw/2+1,NTl);
%     PR = zeros(NC,wndw/2+1,NTr);
%     PC = zeros(NC,wndw/2+1,NTc);
%     for ch = 1:NC
%         ch
%         %Spectrogram computes the short-time fourier transform with data
%         %transformed by the half hamming
%         tic
%         [sS fF tT ~] = spectrogram(Data_R(ch,:),hm,wndw-8,wndw,fs);
%         toc
%         tic
%         %fft data is reorganized like it is online
%         sS2 = [real(single(sS)); flipud(imag(single(sS(1:wndw/2,:))))];
%         clear sS
%         toc
%         %Then the power spectrum is computed by squaring and adding real and
%         %imaginary parts
%
%         normfactor = 1/wndw;
%         mxIndx = wndw/2+1;
%         PS = zeros(wndw/2+1,size(sS2,2));
%         PS(1,:) = sS2(1,:).^2*normfactor;
%         for i = 1:mxIndx-1
%             PS(i+1,:) = (sS2(i+1,:).^2+sS2(wndw+1-i,:).^2)*normfactor;
%         end
%         if mod(wndw,2)==0
%             PS(mxIndx,:)= sS2(mxIndx,:).^2*normfactor;
%         end
%         clear sS2
%         %     PS = log10(PS);
%
%         timeindexL = zeros(size(tT));
%         for i = 1:length(L_trialstart)
%             %Define the time for which the ith left trial is presented.  Cutoff
%             %data a second after the beginning and a second before the end (may
%             %want to fix this later --- define the correct cutoffs in the
%             %definition of L_trialstart/end
%             timeindexL = (tT>(L_trialstart(i)/fs)+1&tT<(L_trialend(i)/fs)-1);
%             PL(ch,:,i) = squeeze(mean(PS(:,timeindexL),2))';
%         end
%         timeindexR = zeros(size(tT));
%         for i = 1:length(R_trialstart)
%             timeindexR = (tT>(R_trialstart(i)/fs)+1&tT<(R_trialend(i)/fs)-1);
%             PR(ch,:,i) = squeeze(mean(PS(:,timeindexR),2))';
%         end
%         timeindexC = zeros(size(tT));
%         for i = 1:length(C_trialstart)
%             timeindexC = (tT>(C_trialstart(i)/fs)+1&tT<(C_trialend(i)/fs)-1);
%             PC(ch,:,i) = squeeze(mean(PS(:,timeindexC),2))';
%         end
%     end
%
%     if Figures_On == 1
%         figure
%         for ch = 1:NC
% %             ttm(2*ch-1) = subplot(4,8,2*ch-1); plot(fF,(squeeze(PL(ch,:,:))),'b'); hold on;
% %             plot(fF,(squeeze(PR(ch,:,:))),'r');
% %             if isempty(PC)==0
% %                 plot(fF,(squeeze(PC(ch,:,:))),'g'); xlim([0 40]);
% %             end
%             ttm(ch) = subplot(4,5,ch); plot(fF,(mean(PL(ch,:,:),3)),'b'); hold on;
%             plot(fF,(mean(PR(ch,:,:),3)),'r');
%             if isempty(PC)==0
%                 plot(fF,(mean(PC(ch,:,:),3)),'g'); xlim([0 40]);
%             end
%         end
%
%
%     end
%
%
%
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     Pdiff = mean(PR-PL,3)./std(PR-PL,[],3);
% %     if size(Data,1)==16
% %         if RerefVal == 2 %Laplacian
% %             Pdiff([1 3],:) = -Pdiff([1 3],:);
% %         elseif RerefVal == 0 || RerefVal == 1 %No reref or CAR
% %             Pdiff = [-Pdiff(4,:); Pdiff(6,:)];
% %         end
% %     elseif size(Data,1)==22
% %
% %     end
% %         Pdiff = mean(Pdiff,1);
% %         [~, maxdiffind] = max(abs(Pdiff(Frange(1):Frange(2))));
% %         if Figures_On == 1
% %             figure
% %             plot(Pdiff);
% %         end
% %         Fdiff = fF(maxdiffind+Frange(1));
% %         %Filter in this band
% %         band_r = [-2 2];
% %         [bb aa] = butter(5,(Fdiff+band_r)/(fs/2));
% %         Data_C = filtfilt(bb,aa,(RefFilt*ArtWeights2*Data)');
% %
% %         if Figures_On == 1
% %             figure
% %             [lpb lpa] = butter(5,.5/(fs/2),'low');
% %             Data_Clp  = filtfilt(lpb,lpa,(Data_C.^2));
% %             plot(Data_Clp(:,[1 2])); hold on; plot(TC,'r')
% %         end
%
%
% end



%% CSP

if DoCSP == 1
    if BuildC == 1
        disp('Not using filtered data for input to CSP')
        Data_C = (RefFilt*ArtWeights2*Data)';
        switch WhichDiff
            case 1
                t_indices = {L_trialstart L_trialend R_trialstart R_trialend};
            case 2
                t_indices = {L_trialstart L_trialend C_trialstart C_trialend};
            case 3
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %6/26/13 -- switched the order of the C's and R's up --
                %hopefully this has the cursor moving in the right
                %direction.
                t_indices = {C_trialstart C_trialend R_trialstart R_trialend};
        end
        [Wp] = CSP_SMR(Data_C,t_indices,NC, NT);
% disp('Using Regularized DSP (W4)!');
    elseif BuildC == 0
    end
elseif DoCSP == 0
    if BuildC == 1
        Wp = eye(size(RefFilt,1));
    elseif BuildC == 0
    end
end

%% Combine Rereferencing, Artifact, CSP matrices into one Spatial Filter
if BuildC == 1
    
    
    if DoCSP == 1
        keepind = [1 NC];
        NC = length(keepind);
    else
        keepind = 1:NC;
    end
    
    Wp = Wp(keepind,:);
    
    SpatFilt = Wp*RefFilt*ArtWeights2;
    %         %Combine Spatial Filter and Artifact Matrix into Spatial Filter for online
    %         %use.  Here i am first doing artifact rejection, then spatial
    %         %filtering.  There is a slight difference with the order.
    %         SpatFilt = ArtWeights'*SpatFilt;
    
    %Save the Spatial Filter - this file is used for running this code in
    %BuildC = 0 mode, I now save a parameter file with this spatial filter
    %for use online.
    try
    dlmwrite([NameBase '\data\' Name Session '\'...
        Name 'S' Session 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
    catch
      dlmwrite(['P:\ALS Proj Data\' Name Session '\'...
        Name 'S' Session 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')    
    end
    %     DataA = ArtWeights2*Data;
    %     DataAF = RefFilt*ArtWeights2*Data;
    %     DataAFC = Wp*RefFilt*ArtWeights2*Data;
    %     DataAFCT = Wp'*RefFilt*ArtWeights2*Data;
    
    
    
    DD = SpatFilt*Data;
    
    
    
elseif BuildC == 0
    try
        SpatFilt = dlmread([NameBaseUserData 'SpatialFilter.txt'],'\t');
    catch
        SpatFilt = dlmread(['P:\ALS Proj Data\' Name CSession '\'...
            Name 'S' CSession 'R' CRun '_SpatialFilter.txt'],'\t');
    end
    

        
    DD = SpatFilt*Data;
    NC = size(DD,1);
    
    if Figures_On==1
        if Art == 2 && DoCSP == 0
            rdata = RefFilt*Data(1:19,:);
            figure
            for i = 1:NC
                nrows = ceil(NC/4);
                ncols = ceil(NC/nrows);
                mm(i) = subplot(nrows,ncols,i); plot(rdata(i,:)); hold on;
                plot(DD(i,:),'r');
            end
            linkaxes(mm,'x');
        end
    end
    
end




% %% Visualize raw,artifacted,filtered,CSPd data
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Visualize time series and power changes
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %separate into left and right
% Llocs = find(StimulusCode>0 & StimulusType == 1);
% Llocs = Clocs(1:StateDuration:end);
% Nlocs = find(StimulusCode>0 & StimulusType == 0);
% Nlocs = Nlocs(1:StateDuration:end);
% CCdata = zeros(length(range),size(Data,2),length(Clocs));
% NNdata = zeros(length(range),size(Data,2),length(Nlocs));
% CCdataA = zeros(length(range),size(DataA,2),length(Clocs));
% NNdataA = zeros(length(range),size(DataA,2),length(Nlocs));
% CCdataAF = zeros(length(range),size(DataAF,2),length(Clocs));
% NNdataAF = zeros(length(range),size(DataAF,2),length(Nlocs));
% CCdataAFC = zeros(length(range),size(DataAFC,2),length(Clocs));
% NNdataAFC = zeros(length(range),size(DataAFC,2),length(Nlocs));
% CCdataAFCT = zeros(length(range),size(DataAFCT,2),length(Clocs));
% NNdataAFCT = zeros(length(range),size(DataAFCT,2),length(Nlocs));
%
% %Plot the spectra of the raw vs rerefed and artifacted channels
% figure
% kk=1;
% for ich =[EEGloc(1) EEGloc(2) EEGloc(end-1-length(EOGloc)) EEGloc(end-length(EOGloc))]
%     subplot(2,2,kk);
%     [Sp1,ff1]=pwelch(Data(:,ich),2*fs,fs,2*fs,fs); hold on;
%     plot(ff1,10*log10(Sp1));
%     [Sp2,ff2]=pwelch(DD(:,ich),2*fs,fs,2*fs,fs);
%     plot(ff2,10*log10(Sp2),'r');
%     kk=kk+1;
% end
%
% DD = DD(:,EEGloc);
% NC = length(EEGloc);
%
%
% % %Plot the raw, filtered, and filtered+artifacted data
% % figure
% % plot(Data(:,1));  hold on;
% % plot(DataF(:,1),'r');
% % if Art == 1
% %     fff = DataF*ArtWeights;
% %     plot(fff(:,1),'g')
% % end
% % plot(DD(:,1),'k');
% % legend('Raw Data','Filtered','Filtered+Artifacted');
% % title('if the green line is not visible, thats good'); clear fff
%
% % %Plot the F and O raw and artifacted, as well as EOG channels
% % figure
% % tt(1)=subplot(311); plot(Data(:,1)); hold on; plot(DD(:,1),'r');
% % tt(2)=subplot(312); plot(Data(:,end-length(EOGloc))); hold on;
% % plot(DD(:,end-length(EOGloc)),'r');
% % tt(3)=subplot(313); plot(Data(:,EOGloc)); hold on
% % linkaxes(tt,'x');
%



%% EyeTracker
if isfield(b,'EyetrackerLeftEyeGazeX')
    P300_EyeTracking_inMIClassifier(GazeX,GazeY, {[L_trialstart...
        L_trialend], [R_trialstart R_trialend], [C_trialstart C_trialend]})
end



%% Calculate Spectra
%The results of PowSpect match the workings of FFTFilter.cpp in the BCI2000
%distribution.  in ::Process(), put code 
%bciout << "Output(" << i << ",0) = " << Output(i,0) << endl;
%This allows us to see the result of the Spectra Calculation.  in BCI2000,
%the input to the fft builds up over time, while the spectrogram function
%matlab starts with a full window of data, and then shifts.  Therefore, the
%output of PowSpect(1,1,1) should match the (wndw/8)th iteration of
%Output(0,0).


wndw = c.FFTWindowLength.NumericValue*fs;
hm = hamming(wndw*2+1);
hm = hm(wndw+1:end-1);
PowSpect = single(zeros(wndw/2+1,(size(DD,2)-wndw)/8+1,NC));
for ch = 1:NC
    ch
    %Spectrogram computes the short-time fourier transform with data
    %transformed by the half hamming
    [sS fF tT ~] = spectrogram(DD(ch,:),hm,wndw-8,wndw,fs);
    
    %fft data is reorganized like it is online
    sS2 = [real(sS); flipud(imag(sS(1:wndw/2,:)))];
    
    %Then the power spectrum is computed by squaring and adding real and
    %imaginary parts
    normfactor = 1/wndw;
    mxIndx = wndw/2+1;
    PS(1,:) = sS2(1,:).^2*normfactor;
    for i = 1:mxIndx-1
        PS(i+1,:) = (sS2(i+1,:).^2+sS2(wndw+1-i,:).^2)*normfactor;
    end
    if mod(wndw,2)==0
        PS(mxIndx,:)= sS2(mxIndx,:).^2*normfactor;
    end
    PowSpect(:,:,ch) = PS;
end
clear sS sS2
%This is commented because there is currently no implementation of log
%spectra in BCI2000
%PowSpect = log10(PowSpect);


% % Visualization of Filtered Data and Spectra
% timeindexL = zeros(size(tT));
% PL = zeros(NC,size(PowSpect,1),NTl);
% % ZL = zeros(NC,size(Z,1),NT/2);
% for i = 1:length(L_trialstart)
%     %Define the time for which the ith left trial is presented.  Cutoff
%     %data a second after the beginning and a second before the end (may
%     %want to fix this later --- define the correct cutoffs in the
%     %definition of L_trialstart/end.  11/15/13 -- Altered this because the
%     %definition of tT starts at .5 seconds.  Therefore the window with 1
%     %second buffers was actually shifted forward .5 seconds.  Shifted back
%     %.5 seconds by taking .5 seconds after the start to 1.5 seconds before
%     %the end.
%     timeindexL = (tT>(L_trialstart(i)/fs)+.5&tT<(L_trialend(i)/fs)-1.5);
%     PL(:,:,i) = squeeze(mean(PowSpect(:,timeindexL,:),2))';
%     %     ZL(:,:,i) = squeeze(mean(Z(:,timeindexL,:),2))';
% end
% timeindexR = zeros(size(tT));
% PR = zeros(NC,size(PowSpect,1),NTr);
% % ZR = zeros(NC,size(Z,1),NT/2);
% for i = 1:length(R_trialstart)
%     timeindexR = (tT>(R_trialstart(i)/fs)+.5&tT<(R_trialend(i)/fs)-1.5);
%     PR(:,:,i) = squeeze(mean(PowSpect(:,timeindexR,:),2))';
%     %     ZR(:,:,i) = squeeze(mean(Z(:,timeindexR,:),2))';
% end
% 
% timeindexC = zeros(size(tT));
% PC = zeros(NC,size(PowSpect,1),NTc);
% % ZR = zeros(NC,size(Z,1),NT/2);
% for i = 1:length(C_trialstart)
%     timeindexC = (tT>(C_trialstart(i)/fs)+.5&tT<(C_trialend(i)/fs)-1.5);
%     PC(:,:,i) = squeeze(mean(PowSpect(:,timeindexC,:),2))';
%     %     ZR(:,:,i) = squeeze(mean(Z(:,timeindexR,:),2))';
% end
% %
% % for ch = 1:NC
% % %     figure(100)
% % %     tt(ch) = subplot(4,4,ch); imagesc(tT,fF,squeeze(Z(:,:,ch))); hold on;
% % %     plot((1:length(TC))/fs,20*TC,'k','LineWidth',2);
% %     figure(101)
% %     xx(ch) = subplot(4,4,ch); imagesc(tT,fF,log10(squeeze(PowSpect(:,:,ch))));hold on;
% %     plot((1:length(TC))/fs,20*TC,'k','LineWidth',2);
% % end
% 
% if BuildC == 1
%     if Figures_On == 1
%         figure('Name','Plot of all features from all three cue types');
%         for ch = 1:NC
%             
%             nrows = ceil(NC/4);
%             ncols = ceil(NC/nrows);
%             ttm(2*ch-1) = subplot(nrows,ncols*2,2*ch-1); plot(fF,(squeeze(PL(ch,:,:))),'c'); hold on;
%             plot(fF,(squeeze(PR(ch,:,:))),'m');
%             plot(fF,(squeeze(PC(ch,:,:))),'g');
%             xlim([0 40]);
%             
%             
%             ttm(2*ch) = subplot(nrows,ncols*2,2*ch); errorbar(fF,mean(PL(ch,:,:),3),std(PL(ch,:,:),[],3),'c'); hold on;
%             errorbar(fF,mean(PR(ch,:,:),3),std(PR(ch,:,:),[],3),'m');
%             errorbar(fF,mean(PC(ch,:,:),3),std(PC(ch,:,:),[],3),'g');
%             xlim([0 40]);
%             %     plot(fF,(mean(ZL(ch,:,:),3)),'b');
%             %     plot(fF,(mean(ZR(ch,:,:),3)),'r');
%         end
%     end
% end


% %Plot for NIH grant -- ref -- 2, csp on, no art
% %saved as NIHgrantSMR-an-sm.fig
% %S1 - AndrewDash_14001 R02 -- csp chan 1 - using left vs right hand imagery
% %S2 - SumithraDash001 R05 -- csp chan 4 - using left hand vs foot imagery
% if Figures_On == 1
%     ch = 1;
%     figure
%     subplot(122);
%     hold on;
%     plot(fF,(mean(PL(ch,:,:),3)),'b'); hold on;
%     plot(fF,(mean(PR(ch,:,:),3)),'r'); plot(fF,(mean(PC(ch,:,:),3)),'g');
%     xlim([0 40]); xlabel('Frequency (Hz)'); ylabel('Power (\muV^2)');
% end





%% Build Classifier
Ctypes = {'SW','L','SS','LO'};
% Findex = fF>=Frange(1) & fF<=Frange(2);
Findex = ismember(fF ,Frange(1):2:Frange(2));
Pfeat = [];
class =[];
PfeatAvg = [];
classAvg = [];
trlct = [];

% if Cfier == 3
%     rlab = 1; llab = 2;
% else
%     rlab = -1; llab = 1;
% end



if Cfier == 3
    TwoClass = [2 1];
else
    TwoClass = [-1 1];
end

switch WhichDiff
    case 1
        TwoCues = {L_trialstart, L_trialend; R_trialstart R_trialend};
    case 2
        TwoCues = {L_trialstart, L_trialend; C_trialstart, C_trialend};
    case 3
        TwoCues = {C_trialstart, C_trialend; R_trialstart, R_trialend};
end

trlctr = 1;
for cc = 1:size(TwoCues,1)
    for tt = 1:length(TwoCues{cc})
        Pfeatch = [];
        %Use 1 second after trial start to 1 second before trial end.
        %11/15/13 -- Altered this because the
    %definition of tT starts at .5 seconds.  Therefore the window with 1
    %second buffers was actually shifted forward .5 seconds.  Shifted back
    %.5 seconds by taking .5 seconds after the start to 1.5 seconds before
    %the end.
        ff = tT>TwoCues{cc,1}(tt)/fs+.5 & tT<TwoCues{cc,2}(tt)/fs-1.5;
        sum(ff)
        for ch = 1:NC
            Pfeatch = [Pfeatch squeeze(PowSpect(Findex,ff,ch))'];
        end
        Pfeat = [Pfeat; Pfeatch];
        PfeatAvg = [PfeatAvg; mean(Pfeatch,1)];
        
        class = [class; TwoClass(cc)*ones(size(Pfeatch,1),1)];
        classAvg = [classAvg; TwoClass(cc)];
        trlct = [trlct; repmat(trlctr,size(Pfeatch,1),1)];
        trlctr = trlctr+1;
    end
end


%Changed to above as of 6/26/13 ---

% for tr = 1:NTl
%
%         switch WhichDiff
%             case 1
%                 ffl = tT*fs>L_trialstart(tr) & tT*fs<L_trialend(tr);
%                 ffr = tT*fs>R_trialstart(tr) & tT*fs<R_trialend(tr);
%             case 2
%                 ffl = tT*fs>L_trialstart(tr) & tT*fs<L_trialend(tr);
%                 ffr = tT*fs>C_trialstart(tr) & tT*fs<C_trialend(tr);
%             case 3
%                 ffl = tT*fs>R_trialstart(tr) & tT*fs<R_trialend(tr);
%                 ffr = tT*fs>C_trialstart(tr) & tT*fs<C_trialend(tr);
%         end
%
%     Pfeatchl = [];
%     Pfeatchr = [];
%     for ch = 1:NC
%         Pfeatchl = [Pfeatchl squeeze(PowSpect(Findex,ffl,ch))'];
%         Pfeatchr = [Pfeatchr squeeze(PowSpect(Findex,ffr,ch))'];
%     end
%     Pfeat = [Pfeat; Pfeatchl; Pfeatchr];
%     PfeatAvg = [PfeatAvg; mean(Pfeatchl,1); mean(Pfeatchr,1)];
%
%     class = [class; llab*ones(size(Pfeatchl,1),1);...
%         rlab*ones(size(Pfeatchr,1),1)];
%     classAvg = [classAvg; llab; rlab];
%     %         class = [class; llab/sum(ffl)*ones(size(Pfeatchl,1),1);...
%     %             rlab/sum(ffl)*ones(size(Pfeatchr,1),1)];
%     %         classAvg = [classAvg; llab/sum(ffl); rlab/sum(ffl)];
%
%     trlct = [trlct; repmat(2*(tr-1)+1,size(Pfeatchl,1),1); ...
%         repmat(2*tr,size(Pfeatchr,1),1)];
%
% end

if Figures_On==1
    figure('Name','Classifier features');
    if AvgF == 0
        for i = 1:NC
            nrows = ceil(NC/4);
            ncols = ceil(NC/nrows);
            
            subplot(nrows,2*ncols,2*i-1); plot(fF(Findex),Pfeat(class==TwoClass(2),...
                (i-1)*sum(Findex)+1:i*sum(Findex)),'r');
            hold on;
            plot(fF(Findex),Pfeat(class==TwoClass(1),...
                (i-1)*sum(Findex)+1:i*sum(Findex)),'b');
            
            subplot(nrows,2*ncols,2*i); errorbar(fF(Findex),mean(Pfeat(class==TwoClass(2),...
                (i-1)*sum(Findex)+1:i*sum(Findex))),...
                std(Pfeat(class==TwoClass(2), (i-1)*sum(Findex)+1:i*sum(Findex))),'r');
            hold on
            errorbar(fF(Findex),mean(Pfeat(class==TwoClass(1),...
                (i-1)*sum(Findex)+1:i*sum(Findex))),...
                std(Pfeat(class==TwoClass(1), (i-1)*sum(Findex)+1:i*sum(Findex))),'b');
        end
    else
        for i = 1:NC
            nrows = ceil(NC/4);
            ncols = ceil(NC/nrows);
            
            subplot(nrows,2*ncols,2*i-1); plot(fF(Findex),PfeatAvg(classAvg==TwoClass(2),...
                (i-1)*sum(Findex)+1:i*sum(Findex)),'r');
            hold on;
            plot(fF(Findex),PfeatAvg(classAvg==TwoClass(1),...
                (i-1)*sum(Findex)+1:i*sum(Findex)),'b');
            
            subplot(nrows,2*ncols,2*i); errorbar(fF(Findex),mean(PfeatAvg(classAvg==TwoClass(2),...
                (i-1)*sum(Findex)+1:i*sum(Findex))),...
                std(PfeatAvg(classAvg==TwoClass(2), (i-1)*sum(Findex)+1:i*sum(Findex))),'r');
            hold on
            errorbar(fF(Findex),mean(PfeatAvg(classAvg==TwoClass(1),...
                (i-1)*sum(Findex)+1:i*sum(Findex))),...
                std(PfeatAvg(classAvg==TwoClass(1), (i-1)*sum(Findex)+1:i*sum(Findex))),'b');
        end
    end
end

%Detect high power (artifactual) features

if AvgF == 1
    OB2 = zeros(size(PfeatAvg,1),1);
    [ssa ssi] = sort(sum(PfeatAvg,2))
    perc25 = ssa(fix((1/4)*length(ssi))); perc75 = ssa(fix((3/4)*length(ssi)));
    maxb = perc75+3*(perc75-perc25);
    minb = perc25-3*(perc75-perc25);
    OB2([ssi(ssa<minb) ssi(ssa>maxb)]) = 1;
%     for i = 1:size(PfeatAvg,1)
%         OB2(i) = sum([sum(PfeatAvg(i,:))<minb sum(PfeatAvg(i,:))>maxb])>0;
%     end
else
    OB2 = zeros(size(Pfeat,1),1);
    [ssa ssi] = sort(sum(Pfeat,2));
    %This method removes too many windows
%     tta = repmat(mean(Pfeat,1),size(Pfeat,1),1);
%     ttv = repmat(var(Pfeat,[],1),size(Pfeat,1),1);
%     tmpp = sum(Pfeat<(tta-3*ttv) | Pfeat>(tta+3*ttv),2)>0;
    perc25 = ssa(fix((1/4)*length(ssi))); perc75 = ssa(fix((3/4)*length(ssi)));
    maxb = perc75+3*(perc75-perc25);
    minb = perc25-3*(perc75-perc25);
    OB2([ssi(ssa<minb) ssi(ssa>maxb)]) = 1;
end
disp(['Windowed spectra removed: ' num2str(sum(OB2))]);
% disp(['Alt Windowed spectra removed: ' num2str(sum(tmpp))]);

% if Figures_On==1
if NC<=6
    figure('Name','Classifier features with outliers highlighted');
    if AvgF==0
        for i = 1:NC
            %Made a plotting change to visualize the windowed spectra that
            %are removed.  Did not make this plotting change for the AvgF =
            %1 case.
            nrows = ceil(NC/4);
            ncols = ceil(NC/nrows);
            
            subplot(nrows,2*ncols,2*i-1); 
            temp2 = find(class==TwoClass(2));
            tempA = find(~OB2);
            plot(fF(Findex),Pfeat(temp2(ismember(temp2,tempA)),...
                (i-1)*sum(Findex)+1:i*sum(Findex)),'Color',[1 .9 .9]);
            hold on;
            temp1 = find(class==TwoClass(1));
            plot(fF(Findex),Pfeat(temp1(ismember(temp1,tempA)),...
                (i-1)*sum(Findex)+1:i*sum(Findex)),'Color',[.9 .9 1]);
            
            tempK = find(OB2);
            if sum(ismember(temp2,tempK))~=0
            plot(fF(Findex),Pfeat(temp2(ismember(temp2,tempK)),...
                (i-1)*sum(Findex)+1:i*sum(Findex)),'r');
            end;
            hold on;
            if sum(ismember(temp1,tempK))~=0
            plot(fF(Findex),Pfeat(temp1(ismember(temp1,tempK)),...
                (i-1)*sum(Findex)+1:i*sum(Findex)),'b');
            end;
            
            subplot(nrows,2*ncols,2*i); 
            errorbar(fF(Findex),mean(Pfeat(temp2,...
                (i-1)*sum(Findex)+1:i*sum(Findex))),...
                std(Pfeat(temp2, (i-1)*sum(Findex)+1:i*sum(Findex))),'Color',[1 .9 .9]);
            hold on;
            errorbar(fF(Findex),mean(Pfeat(temp1,...
                (i-1)*sum(Findex)+1:i*sum(Findex))),...
                std(Pfeat(temp1, (i-1)*sum(Findex)+1:i*sum(Findex))),'Color',[.9 .9 1]);
            errorbar(fF(Findex),mean(Pfeat(temp2(ismember(temp2,tempA)),...
                (i-1)*sum(Findex)+1:i*sum(Findex))),...
                std(Pfeat(temp2(ismember(temp2,tempA)), (i-1)*sum(Findex)+1:i*sum(Findex))),'r');
            errorbar(fF(Findex),mean(Pfeat(temp1(ismember(temp1,tempA)),...
                (i-1)*sum(Findex)+1:i*sum(Findex))),...
                std(Pfeat(temp1(ismember(temp1,tempA)), (i-1)*sum(Findex)+1:i*sum(Findex))),'b');
        end
    else
        for i = 1:NC
            nrows = ceil(NC/4);
            ncols = ceil(NC/nrows);
            
            subplot(nrows,2*ncols,2*i-1); plot(fF(Findex),PfeatAvg(classAvg==TwoClass(2),...
                (i-1)*sum(Findex)+1:i*sum(Findex)),'r');
            hold on;
            plot(fF(Findex),PfeatAvg(classAvg==TwoClass(1),...
                (i-1)*sum(Findex)+1:i*sum(Findex)),'b');
            
            subplot(nrows,2*ncols,2*i); errorbar(fF(Findex),mean(PfeatAvg(classAvg==TwoClass(2),...
                (i-1)*sum(Findex)+1:i*sum(Findex))),...
                std(PfeatAvg(classAvg==TwoClass(2), (i-1)*sum(Findex)+1:i*sum(Findex))),'r');
            hold on
            errorbar(fF(Findex),mean(PfeatAvg(classAvg==TwoClass(1),...
                (i-1)*sum(Findex)+1:i*sum(Findex))),...
                std(PfeatAvg(classAvg==TwoClass(1), (i-1)*sum(Findex)+1:i*sum(Findex))),'b');
        end
    end
end
% end


%Remove high power (artifactual) features
if AvgF == 1
    PfeatAvg = PfeatAvg(~OB2,:);
    classAvg = classAvg(~OB2);
else
    Pfeat = Pfeat(~OB2,:);
    class = class(~OB2);
end
trlct = trlct(~OB2);

if BuildC == 1
    
    
    %Open Test Parameter File for writing
    fpm = fopen([NameBase '\parms\ALSDash\FINAL_Dasher_test_22_3EOG_2amp.prm'],'r');
    fpm2=textscan(fpm,'%s','delimiter','\n');
    fclose(fpm);
    %Open Free parameter File for writing
    fem = fopen([NameBase '\parms\ALSDash\FINAL_Dasher_free_22_3EOG_2amp.prm'],'r');
    fem2=textscan(fem,'%s','delimiter','\n');
    fclose(fem);
    for ll = 1:length(fpm2{1})
        LinInd(ll) = ~isempty(strfind(fpm2{1}{ll},'Filtering:Linear'));
        SNInd(ll) = ~isempty(strfind(fpm2{1}{ll},'SubjectName'));
        SEInd(ll) = ~isempty(strfind(fpm2{1}{ll},'SubjectSession'));
        SFInd(ll) = ~isempty(strfind(fpm2{1}{ll},'SpatialFilter matrix'));
        FIInd(ll) = ~isempty(strfind(fpm2{1}{ll},'FFTInputChannels'));
        NOInd(ll) = ~isempty(strfind(fpm2{1}{ll},'NormalizerOffsets'));
        NGInd(ll) = ~isempty(strfind(fpm2{1}{ll},'NormalizerGains'));
    end
    LinInd = find(LinInd==1);
    SNInd = find(SNInd==1);
    SEInd = find(SEInd==1);
    SFInd = find(SFInd==1);
    FIInd = find(FIInd==1);
    NOInd = find(NOInd==1);
    NGInd = find(NGInd==1);
    SBegInd = strfind(fpm2{1}{SNInd},'='); SEndInd = strfind(fpm2{1}{SNInd},'Name %');
    SEBegInd = strfind(fpm2{1}{SEInd},'='); SEEndInd = strfind(fpm2{1}{SEInd},'% %');
    SFBegInd = strfind(fpm2{1}{SFInd},'='); SFEndInd = strfind(fpm2{1}{SFInd},'//');
    BegInd = strfind(fpm2{1}{LinInd},'}'); EndInd = strfind(fpm2{1}{LinInd},'//');
    BegInd2 = strfind(fpm2{1}{LinInd},'= '); EndInd2 = strfind(fpm2{1}{LinInd},' {');
    FIBegInd = strfind(fpm2{1}{FIInd},'='); FIEndInd = strfind(fpm2{1}{FIInd},'//');
    NOBegInd = strfind(fpm2{1}{NOInd},'='); NOEndInd = strfind(fpm2{1}{NOInd},'% %');
    NGBegInd = strfind(fpm2{1}{NGInd},'='); NGEndInd = strfind(fpm2{1}{NGInd},'% %');
    
    for ll = 1:length(fem2{1})
        LinIndf(ll) = ~isempty(strfind(fem2{1}{ll},'Filtering:Linear'));
        SNIndf(ll) = ~isempty(strfind(fem2{1}{ll},'SubjectName'));
        SEIndf(ll) = ~isempty(strfind(fem2{1}{ll},'SubjectSession'));
        SFIndf(ll) = ~isempty(strfind(fem2{1}{ll},'SpatialFilter matrix'));
        FIIndf(ll) = ~isempty(strfind(fem2{1}{ll},'FFTInputChannels'));
        NOIndf(ll) = ~isempty(strfind(fem2{1}{ll},'NormalizerOffsets'));
        NGIndf(ll) = ~isempty(strfind(fem2{1}{ll},'NormalizerGains'));
    end
    LinIndf = find(LinIndf==1);
    SNIndf = find(SNIndf==1);
    SEIndf = find(SEIndf==1);
    SFIndf = find(SFIndf==1);
    FIIndf = find(FIIndf==1);
    NOIndf = find(NOIndf==1);
    NGIndf = find(NGIndf==1);
    SBegIndf = strfind(fem2{1}{SNIndf},'='); SEndIndf = strfind(fem2{1}{SNIndf},'Name %');
    SEBegIndf = strfind(fem2{1}{SEIndf},'='); SEEndIndf = strfind(fem2{1}{SEIndf},'% %');
    SFBegIndf = strfind(fem2{1}{SFIndf},'='); SFEndIndf = strfind(fem2{1}{SFIndf},'//');
    BegIndf = strfind(fem2{1}{LinIndf},'}'); EndIndf = strfind(fem2{1}{LinIndf},'//');
    BegInd2f = strfind(fem2{1}{LinIndf},'= '); EndInd2f = strfind(fem2{1}{LinIndf},' {');
    FIBegIndf = strfind(fem2{1}{FIIndf},'='); FIEndIndf = strfind(fem2{1}{FIIndf},'//');
    NOBegIndf = strfind(fem2{1}{NOIndf},'='); NOEndIndf = strfind(fem2{1}{NOIndf},'% %');
    NGBegIndf = strfind(fem2{1}{NGIndf},'='); NGEndIndf = strfind(fem2{1}{NGIndf},'% %');
    
    
    if Cfier == 1
        if AvgF == 0
            [B, SE, PVAL, INMODEL, STATS] = stepwisefit(Pfeat,class,...
                'penter',.001,'premove',.01);
        elseif AvgF == 1
            [B, SE, PVAL, INMODEL, STATS] = stepwisefit(PfeatAvg,classAvg,...
                'penter',.001,'premove',.01);
        end
        Classweight = -B(INMODEL==1);
        SigLocs = find(INMODEL==1);
        
        %         [~, keepind] = sort(PVAL);
        %         if length(keepind)>10
        %             keepind = keepind(1:10);
        %         end
        %         SigLocs = keepind;
        %         Pfeat2 = Pfeat(:,SigLocs);
        %         [B, SE, PVAL, INMODEL, STATS] = stepwisefit(Pfeat2,class);
        %         Classweight = -B(INMODEL==1);
        %         SigLocs = SigLocs(INMODEL==1);
        
        
    elseif Cfier == 2
        if AvgF == 0
            [clsv,err,post,logp,coeff] = classify(Pfeat,Pfeat,class);
        elseif Avv == 1
            [clsv,err,post,logp,coeff] = classify(PfeatAvg,PfeatAvg,classAvg);
        end
        Classweight = coeff(1,2).linear;
        SigLocs = find(Classweight~=0);
        
    elseif Cfier == 3
        if AvgF == 0
            cd (['C:\Documents and Settings\' pcusr '\Google Drive\Year 4\'...
                'ALS Project\Part2 - P300\P300Classifier\LDA Code\Redo'])
            [Gam,z,SigLocs] = runLDA2([],Pfeat,class,10);
            cd (['C:\Documents and Settings\' pcusr '\Google Drive\Year 4\'...
                'ALS Project\Part2 - P300\SMRClassifier'])
        elseif AvgF == 1
            cd (['C:\Documents and Settings\' pcusr '\Google Drive\Year 4\'...
                'ALS Project\Part2 - P300\P300Classifier\LDA Code\Redo'])
            [Gam,z,SigLocs] = runLDA2([],PfeatAvg,classAvg,10);
            cd (['C:\Documents and Settings\' pcusr '\Google Drive\Year 4\'...
                'ALS Project\Part2 - P300\SMRClassifier'])
        end
        
        %Sometimes the polarity is reversed on the output.  We always want
        %left feedback to be negative.  In runlda2 -- the groups are
        %organized so that the first half of z has a label of 1
        %(right), and the second half of z has a label of 2 (left).
        if mean(z(1:length(z)/2))-mean(z(length(z)/2+1:end))>0
            Classweight = Gam;
        else
            Classweight = -Gam;
        end
        
    elseif Cfier == 4 %LASSO
        if AvgF == 0
            [B2, STATS2] =  lasso(Pfeat,class,'CV',10);
        elseif AvgF == 1
            [B2, STATS2] =  lasso(PfeatAvg,classAvg,'CV',10);
        end
        %Retain variables in the model corresponding to the sparsest
        %model within one standard error of the minimum cross validated
        %MSE
        Classweight = -B2(:,STATS2.Index1SE);
        SigLocs = find(Classweight~=0);
        Classweight = Classweight(SigLocs);
        
    else
        disp('Wrong Classifier')
        return
    end
    
    
    CM = zeros(length(SigLocs),4);
    %     CM = repmat(cellstr(''), length(SigLocs),4);
    if isempty(CM)
        disp('No significant features!');
        return
    else
        %         CM(:,1) = mat2cell(ceil(SigLocs./sum(Findex)),1,ones(size(SigLocs,2),1))';
        CM(:,1) = (ceil(SigLocs./sum(Findex))');
        Ctms = mod(SigLocs,sum(Findex));  Ctms(Ctms==0)=sum(Findex);
        Ffd = find(Findex==1);
        %         CM(:,2) = mat2cell(fF(Ffd(Ctms)),ones(size(SigLocs,2),1),1); %Input element (sample #)
        CM(:,2) = fF(Ffd(Ctms)); %Input element (sample #)
        %         Frange(1)-1+SigLocs-(length(Frange).*(CM(:,1)-1))';
        CM(:,3) = 1;  %Moves cursor left and right
        CM(:,4) = Classweight;
    end
    
    %Looking at the expected classification for trial averages
    ffrt = fF(Findex);
    for ii = 1:size(PfeatAvg,1)
        cursA(ii) = 0;
        for jj = 1:size(CM,1)
            chnn = CM(jj,1);
            ffr = (chnn-1)*sum(Findex)+1:chnn*sum(Findex);
            ffrx = ffrt==CM(jj,2);
            cursA(ii) = cursA(ii) + CM(jj,4)*PfeatAvg(ii,ffr(ffrx));
        end
    end
    
    %Looking at the expected classification for individual segments
    for ii = 1:size(Pfeat,1)
        curs(ii) = 0;
        for jj = 1:size(CM,1)
            chnn = CM(jj,1);
            chind = (chnn-1)*length(ffrt)+1:chnn*length(ffrt);
            ffrx = ffrt==CM(jj,2);
            curs(ii) = curs(ii) + CM(jj,4)*Pfeat(ii,chind(ffrx));
        end
    end
    for tt = unique(trlct)'
        trlindex = trlct==tt;
        PfeatAvg(tt,:)=mean(Pfeat(trlindex,:));
        cursA2(tt) = mean(curs(trlindex));
    end
     

    %What is optimal decision boundary?
    %Set the gain and the offset of the Normalizer
    
    lcurs = curs(class==TwoClass(1));
    rcurs = curs(class==TwoClass(2));
    k=1;
    Offrange = min(rcurs):.01:max(lcurs);
    if isempty(Offrange)
        NormalizerOffset = (min(rcurs)+max(lcurs))/2;
    else
        for thr = Offrange %Cycle through possible thresholds
            psums(k) = sum([lcurs<thr rcurs>thr])/(length(lcurs)+length(rcurs));
            k=k+1;
        end
        [ma mi] = max(psums);
        NormalizerOffset = Offrange(mi);
    end
    NormalizerGain = 1/std(curs);
    
%     trlindex = trlct==1;
%     ccheck = curs(trlindex)-NormalizerOffset;
%     ccheck = cumsum(ccheck);
%     figure
%     plot(ccheck);
%     
    
    
    if Figures_On == 1
        figure
        subplot(211);
        hist(cursA(classAvg==TwoClass(1)),10);
        hold on;
        hist(cursA(classAvg==TwoClass(2)),10);
        h = findobj(gca,'Type','patch');
        set(h(1),'FaceColor','r','EdgeColor','w','FaceAlpha',.5);
        %weird that this labeling of h(1) and h(2) is backwards
        set(h(2),'FaceColor','b','EdgeColor','w','FaceAlpha',.5);
        subplot(212);
        hist(curs(class==TwoClass(1)),100);
        hold on;
        hist(curs(class==TwoClass(2)),100);
        h = findobj(gca,'Type','patch');
        set(h(1),'FaceColor','r','EdgeColor','w','FaceAlpha',.5);
        set(h(2),'FaceColor','b','EdgeColor','w','FaceAlpha',.5);
        line([NormalizerOffset NormalizerOffset],[0 100],'Color','k','LineWidth',3)
    end
    
    
    %Instead of the line on the histogram plot, do an ROC curve
    %If we consider "truth" as being a left trial
    
    if isempty(Offrange)
    else
        k=1;
        for thr = Offrange
            tp(k) = sum(lcurs<thr); fn(k) = sum(lcurs>thr);
            tn(k) = sum(rcurs>thr); fp(k) = sum(rcurs<thr); k = k+1;
        end
        tpr = tp./(tp+fn);
        fpr = fp./(fp+tn);
        auckinda = tpr.*(1-fpr);
        [~, emaxi] = max(auckinda);
        Offrange(emaxi);
        lineslope = 1;
        range2 = 0:.01:1;
        deg45line = lineslope*(range2-fpr(emaxi))+tpr(emaxi);
        if Figures_On == 1
            figure
            plot(fpr,tpr); hold on
            plot(range2,deg45line,'r'); xlim([0 1]); ylim([0 1]);
        end
        
        %Do the two methods produce the same result?
        if mi==emaxi
            disp('Yay!');
        else
            disp('Different :(');
        end
    end
    
    
    
    
    %Change the classifier matrix to reflect this boundary
    
    
    
    
    
    
    %%Append Hz and convert to savable form
    CM2 = mat2cell(CM,ones(size(CM,1),1),ones(size(CM,2),1));
    CM2(:,2) = cellfun(@num2str,CM2(:,2),'UniformOutput',false);
    for ii = 1:size(CM2,1);
        CM2{ii,2} = [CM2{ii,2} 'Hz'];
    end
    %     fid=fopen(['C:\Documents and Settings\amg5106\My Documents\'...
    %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' Name Session '\'...
    %         Name 'S' Session 'R' NewRunName '_' Ctypes{Cfier} 'Classifier.txt'],'wt');
    fid=fopen([NameBase '\data\' Name Session '\'...
        Name 'S' Session 'R' NewRunName '_' Ctypes{Cfier} 'Classifier.txt'],'wt');
    if fid==-1
        fid=fopen(['P:\\ALS Proj Data\' Name Session '\'...
            Name 'S' Session 'R' NewRunName '_' Ctypes{Cfier} 'Classifier.txt'],'wt');
    end
    [rows,cols]=size(CM2)
    for i=1:rows
        fprintf(fid,'%d\t',CM2{i,1});
        fprintf(fid,'%s\t',CM2{i,2});
        fprintf(fid,'%d\t',CM2{i,3});
        fprintf(fid,'%e\n',CM2{i,4});
    end
    fprintf(fid,'%f\t%d\t%d\t%d\n',NormalizerGain,0,0,0);
    fprintf(fid,'%f\t%d\t%d\t%d\n',NormalizerOffset,0,0,0);
    fclose(fid);
    
    
    %Write to parameter files
    
    %Save Test parameter file
    %Change Subject Name
    fpm2{1}{SNInd} = [fpm2{1}{SNInd}(1:SBegInd+1) Name...
        fpm2{1}{SNInd}(SEndInd-1:end)];
    %Change Session
    fpm2{1}{SEInd} = [fpm2{1}{SEInd}(1:SEBegInd+1) Session...
        fpm2{1}{SEInd}(SEEndInd-1:end)];
    %Change Spatial Filter
    prmsf = [num2str(size(SpatFilt))  ' ' num2str(reshape(SpatFilt',1,...
        size(SpatFilt,1)*size(SpatFilt,2)))];
    fpm2{1}{SFInd} = [fpm2{1}{SFInd}(1:SFBegInd+1) prmsf...
        fpm2{1}{SFInd}(SFEndInd-1:end)];
    %Change Classifier
    CM2 = cellfun(@(x) [num2str(x) ' '],CM2,'UniformOutput',false);
    prmcls = reshape(CM2',1,size(CM2,1)*size(CM2,2));
    fpm2{1}{LinInd} = [fpm2{1}{LinInd}(1:BegInd2+1) num2str(size(CM,1)) ...
        fpm2{1}{LinInd}(EndInd2:BegInd+1) prmcls...
        fpm2{1}{LinInd}(EndInd-1:end)];
    fpm2{1}{LinInd} = horzcat(fpm2{1}{LinInd}{:});
    %Change FFT Input Channels
    fpm2{1}{FIInd} = [fpm2{1}{FIInd}(1:FIBegInd+1) num2str(size(SpatFilt,1)) ...
        ' ' num2str(1:size(SpatFilt,1)) fpm2{1}{FIInd}(FIEndInd-1:end)];
    %Change Normalizer Offsets
    fpm2{1}{NOInd} = [fpm2{1}{NOInd}(1:NOBegInd+1) ...
        num2str([2 NormalizerOffset 0 0]) fpm2{1}{NOInd}(NOEndInd-1:end)];
    %Change Normalizer Gains
    fpm2{1}{NGInd} = [fpm2{1}{NGInd}(1:NGBegInd+1) ...
        num2str([2 NormalizerGain 1 0]) fpm2{1}{NGInd}(NGEndInd-1:end)];
    
    try
        dlmcell([NameBase '\data\' Name Session '\TEST_ParamFile'...
            Name 'S' Session 'R' NewRunName '_' Ctypes{Cfier} '.prm'],fpm2{1});
    catch
        dlmcell(['P:\\ALS Proj Data\' Name Session '\TEST_ParamFile'...
            Name 'S' Session 'R' NewRunName '_' Ctypes{Cfier} '.prm'],fpm2{1});
    end
    
    %Save Free parameter file
    %Change Subject Name
    fem2{1}{SNIndf} = [fem2{1}{SNIndf}(1:SBegIndf+1) Name...
        fem2{1}{SNIndf}(SEndIndf-1:end)];
    %Change Session
    fem2{1}{SEIndf} = [fem2{1}{SEIndf}(1:SEBegIndf+1) Session...
        fem2{1}{SEIndf}(SEEndIndf-1:end)];
    %Change Spatial Filter
    prmsf = [num2str(size(SpatFilt))  ' ' num2str(reshape(SpatFilt',1,...
        size(SpatFilt,1)*size(SpatFilt,2)))];
    fem2{1}{SFIndf} = [fem2{1}{SFIndf}(1:SFBegIndf+1) prmsf...
        fem2{1}{SFIndf}(SFEndIndf-1:end)];
    %Change Classifier
    CM2 = cellfun(@(x) [num2str(x) ' '],CM2,'UniformOutput',false);
    prmcls = reshape(CM2',1,size(CM2,1)*size(CM2,2));
    fem2{1}{LinIndf} = [fem2{1}{LinIndf}(1:BegInd2f+1) num2str(size(CM,1)) ...
        fem2{1}{LinIndf}(EndInd2f:BegIndf+1) prmcls...
        fem2{1}{LinIndf}(EndIndf-1:end)];
    fem2{1}{LinIndf} = horzcat(fem2{1}{LinIndf}{:});
    %Change FFT Input Channels
    fem2{1}{FIIndf} = [fem2{1}{FIIndf}(1:FIBegIndf+1) num2str(size(SpatFilt,1)) ...
        ' ' num2str(1:size(SpatFilt,1)) fem2{1}{FIIndf}(FIEndIndf-1:end)];
    %Change Normalizer Offsets
    fem2{1}{NOIndf} = [fem2{1}{NOIndf}(1:NOBegIndf+1) ...
        num2str([2 NormalizerOffset 0 0]) fem2{1}{NOIndf}(NOEndIndf-1:end)];
    %Change Normalizer Gains
    fem2{1}{NGIndf} = [fem2{1}{NGIndf}(1:NGBegIndf+1) ...
        num2str([2 NormalizerGain 1 0]) fem2{1}{NGIndf}(NGEndIndf-1:end)];
    
    try
    dlmcell([NameBase '\data\' Name Session '\FREE_ParamFile'...
        Name 'S' Session 'R' NewRunName '_' Ctypes{Cfier} '.prm'],fem2{1});
    catch
        dlmcell(['P:\\ALS Proj Data\' Name Session '\FREE_ParamFile'...
            Name 'S' Session 'R' NewRunName '_' Ctypes{Cfier} '.prm'],fem2{1});
    end
    
    
elseif BuildC == 0 %Load other classifier Matrix
    
    try
        [Cinp Cfreq Cout Cwei] = textread([NameBaseUserData Ctypes{Cfier} 'Classifier.txt'],'%f %f %*s %f %f')
    catch
        [Cinp Cfreq Cout Cwei] = textread(['P:\\ALS Proj Data\' Name CSession '\'...
            Name 'S' CSession 'R' CRun '_' Ctypes{Cfier} 'Classifier.txt'],'%f %f %*s %f %f');
    end
    
    NormalizerGain = Cinp(end-1,1);
    NormalizerOffset = Cinp(end,1);
    CM = [Cinp Cfreq Cout Cwei];
    CM = CM(1:end-2,:)
    
    %Looking at the expected classification for trial averages
    ffrt = fF(Findex);
    for ii = 1:size(PfeatAvg,1)
        cursA(ii) = 0;
        for jj = 1:size(CM,1)
            chnn = CM(jj,1);
            ffr = (chnn-1)*sum(Findex)+1:chnn*sum(Findex);
            ffrx = ffrt==CM(jj,2);
            cursA(ii) = cursA(ii) + CM(jj,4)*PfeatAvg(ii,ffr(ffrx));
        end
    end
    
    %Looking at the expected classification for individual segments
    for ii = 1:size(Pfeat,1)
        curs(ii) = 0;
        for jj = 1:size(CM,1)
            chnn = CM(jj,1);
            chind = (chnn-1)*length(ffrt)+1:chnn*length(ffrt);
            ffrx = ffrt==CM(jj,2);
            curs(ii) = curs(ii) + CM(jj,4)*Pfeat(ii,chind(ffrx));
        end
        curs(ii) = curs(ii);
    end
    for tt = unique(trlct)'
        trlindex = trlct==tt;
        PfeatAvg(tt,:)=mean(Pfeat(trlindex,:));
        cursA2(tt) = mean(curs(trlindex));
    end
    
    
%     trlindex = trlct==18;
%     ccheck = curs(trlindex)-NormalizerOffset;
%     ccheck = cumsum(ccheck);
%     figure
%     plot(ccheck);

%     
    %What is optimal decision boundary?
    lcurs = curs(class==TwoClass(1));
    rcurs = curs(class==TwoClass(2));
    k=1;
    Offrange = min(rcurs):.01:max(lcurs);
    %     for thr = Offrange %Cycle through possible thresholds
    %         psums(k) = sum([lcurs<thr rcurs>thr])/(length(lcurs)+length(rcurs));
    %         k=k+1;
    %     end
    %     [ma mi] = max(psums)
    %Set the gain and the offset of the Normalizer
    %     NormalizerOffset = Offrange(mi);
    %     NormalizerGain = 1/var(curs);
    
    
    if Figures_On == 1
        figure
        subplot(211);
        hist(cursA(classAvg==TwoClass(1)),10);
        hold on;
        hist(cursA(classAvg==TwoClass(2)),10);
        h = findobj(gca,'Type','patch');
        set(h(1),'FaceColor','r','EdgeColor','w','FaceAlpha',.5);
        %weird that this labeling of h(1) and h(2) is backwards
        set(h(2),'FaceColor','b','EdgeColor','w','FaceAlpha',.5);
        subplot(212);
        hist(curs(class==TwoClass(1)),100);
        hold on;
        hist(curs(class==TwoClass(2)),100);
        h = findobj(gca,'Type','patch');
        set(h(1),'FaceColor','r','EdgeColor','w','FaceAlpha',.5);
        set(h(2),'FaceColor','b','EdgeColor','w','FaceAlpha',.5);
        line([NormalizerOffset NormalizerOffset],[0 100],'Color','k','LineWidth',3)
    end
    
    
    %     %Instead of the line on the histogram plot, do an ROC curve
    %     %If we consider "truth" as being a left trial
    %     k=1;
    %     for thr = Offrange
    %         tp(k) = sum(lcurs<thr); fn(k) = sum(lcurs>thr);
    %         tn(k) = sum(rcurs>thr); fp(k) = sum(rcurs<thr); k = k+1;
    %     end
    %     tpr = tp./(tp+fn);
    %     fpr = fp./(fp+tn);
    %     auckinda = tpr.*(1-fpr);
    %     [~, emaxi] = max(auckinda);
    %     Offrange(emaxi);
    %     lineslope = 1;
    %     range2 = 0:.01:1;
    %     deg45line = lineslope*(range2-fpr(emaxi))+tpr(emaxi);
    %     if Figures_On == 1
    %         figure
    %         plot(fpr,tpr); hold on
    %         plot(range2,deg45line,'r'); xlim([0 1]); ylim([0 1]);
    %     end
    %
    %     %Do the two methods produce the same result? -- i dont think they
    %     %should... one maximizes the  accuracy (TP+TN)/(TP+TN+FP+FN), and one
    %     %maximizes AUC (kinda).
    %     if mi==emaxi
    %         disp('Yay!');
    %     else
    %         disp('Different :(');
    %     end
else
end

%% New Feedback with new classifier
%Cursor speed (DasherSpeller.cpp) is based on the Feedback interval length
if Figures_On == 1
    figure
    title('Right is 1, left is 2');
    plot(1000*Feedback,'r');hold on;
    plot(999*TC,'b');
    plot(Cursor,'k');
    if BuildC == 0
        Classifiers = {ClassifierUsed, CM};
        Normvalues = {[c.NormalizerGains.NumericValue(1) c.NormalizerOffsets.NumericValue(1)]...
            [NormalizerGain NormalizerOffset]};
    elseif BuildC == 1
        Classifiers = {CM};
        Normvalues = {[NormalizerGain NormalizerOffset]}
    end;
    
    colors = {'g','m'};
    for cv = 1:length(Classifiers)
        disp(['Classifier' num2str(cv)]);
        clsfr = Classifiers{cv};
        NormV = Normvalues{cv}
        CursorSpeed = 100 / (c.FeedbackDuration.NumericValue*32 * 2);
        StartingPos = double([50;50;50]);
        TargetPos = c.Targets.Value(:,1:3);
        LTv = str2double(TargetPos{2,1}); RTv = str2double(TargetPos{1,1});
        CursorO = double(StartingPos(1)*ones(3,size(PowSpect,2)));
        tindex = wndw;
        targethit = 0;
        % for tt =wndw/8+1:95
        for tt = wndw/8+1:size(PowSpect,2)
            %Is feedback running
            if Feedback(tindex) == 1
                %Did feedback just began, reset Cursor to Starting Position
                if Feedback(tindex-8) == 0
                    CursorO(:,tt) = StartingPos;
                    targethit = 0;
                else
                    %Increment Cursor
                    if targethit == 0
                        cursinc= double(zeros(3,1));
                        cursss(tt) = 0;
                        for i = 1:size(clsfr,1)
                            controlsignal = double(clsfr(i,4)*...
                                PowSpect(fF==clsfr(i,2),tt-wndw/8,clsfr(i,1)));
                            cursinc(clsfr(i,3)) = cursinc(clsfr(i,3))...
                                + controlsignal;
                            %                         cursss(tt) = cursss(tt) + clsfr(i,4)*...
                            %                             PowSpect(find(fF==clsfr(i,2)),tt-wndw/8,clsfr(i,1));
                        end
                        %Normalizer
                        cursinc = NormV(1)*(cursinc-NormV(2));
                        
                        CursorO(:,tt) = CursorO(:,tt-1)+CursorSpeed*cursinc;
                        
                 
                        
                        %Check if a target has been hit
                        if WhichDiff == 1 || WhichDiff == 2
                            
                            if TC(tindex) == 1 %right target is displayed
                                if CursorO(1,tt) >= RTv
                                    disp([num2str(tindex) 'in R_t'])
                                    targethit = 1;
                                end
                            elseif TC(tindex) == 2 %left/cross target
                                if CursorO(1,tt) <= LTv
                                    disp([num2str(tindex) 'in L_t'])
                                    targethit = 1;
                                end
                            else
                            end
                        else
                            if TC(tindex) == 1 %right target is displayed
                                if CursorO(1,tt) >= RTv
                                    disp('in tt')
                                    targethit = 1;
                                end
                            elseif TC(tindex) == 2 %cross target
                                if CursorO(1,tt) <= LTv
                                    disp('in bt')
                                    targethit = 1;
                                end
                            else
                            end
                        end

                        %Keep target within bounds
                        if CursorO(1,tt) > 95
                            CursorO(1,tt) = 95;
                        elseif CursorO(1,tt) < 5
                            CursorO(1,tt) = 5;
                        else
                        end
                    else
                        CursorO(:,tt) = CursorO(:,tt-1);
                    end
                end
            else
                CursorO(:,tt) = CursorO(:,tt-1);
            end
            tindex = tindex+8;           
            
        end
        
        coordToState = 40.95;  %This is from the DoFeedback func. of DasherSpeller.cpp
        CursorO = fix(CursorO.*coordToState);
        %Only interested in the second index - the Y position of the cursor
        CursorO = (CursorO(1,:));
        CursorO = repmat(CursorO,8,1);
        Output4 = reshape(CursorO,1,size(CursorO,1)*size(CursorO,2));
        
        
        plot(Output4,'Color',colors{cv})
    end
    
    
    
    %     %% Before Changing pre3/25/13 -- reproduced output
    %     figure
    %     plot(1000*Feedback,'r');hold on;
    %     plot(999*TC,'b');
    %     plot(Cursor,'k');
    %     Classifiers = {ClassifierUsed, CM};
    %     colors = {'g','m'};
    %     for cv = 1:2
    %         disp(['Classifier' num2str(cv)]);
    %         clsfr = Classifiers{cv};
    %         CursorSpeed = 100 / (c.FeedbackDuration.NumericValue*32 * 2);
    %         StartingPos = double([50;50;50]);
    %         TargetPos = [50 50;10 90;50 50];
    %         CursorO = double(zeros(3,size(PowSpect,2)));
    %         tindex = wndw;
    %         targethit = 0;
    %         % for tt =wndw/8+1:95
    %         for tt = wndw/8+1:size(PowSpect,2)
    %             %Is feedback running
    %             if Feedback(tindex) == 1
    %                 %Did feedback just began, reset Cursor to Starting Position
    %                 if Feedback(tindex-8) == 0
    %                     CursorO(:,tt) = StartingPos;
    %                     targethit = 0;
    %                 else
    %                     %Increment Cursor
    %                     if targethit == 0
    %                         cursinc = double(zeros(3,1));
    %                         cursss(tt) = 0;
    %                         for i = 1:size(clsfr,1)
    %                             cursinc(clsfr(i,3)) = cursinc(clsfr(i,3))...
    %                                 + double(CursorSpeed*clsfr(i,4)*...
    %                                 PowSpect(find(fF==clsfr(i,2)),tt-wndw/8,clsfr(i,1)));
    %                             cursss(tt) = cursss(tt) + clsfr(i,4)*...
    %                                 PowSpect(find(fF==clsfr(i,2)),tt-wndw/8,clsfr(i,1));
    %                         end
    %                         CursorO(:,tt) = CursorO(:,tt-1)+cursinc;
    %
    %                         %Check if a target has been hit
    %                         if TC(tindex) == 1 %Top target is displayed
    %                             if CursorO(1,tt) >= TargetPos(2,2)
    %                                 disp('in tt');
    %                                 targethit = 1;
    %                             end
    %                         elseif TC(tindex) == 2
    %                             if CursorO(1,tt) <= TargetPos(2,1)
    %                                 disp('in bt')
    %                                 targethit = 1;
    %                             end
    %                         else
    %                         end
    %
    %                         %Keep target within bounds
    %                         if CursorO(1,tt) > 95
    %                             CursorO(1,tt) = 95;
    %                         elseif CursorO(1,tt) < 5
    %                             CursorO(1,tt) = 5;
    %                         else
    %                         end
    %                     else
    %                         CursorO(:,tt) = CursorO(:,tt-1);
    %                     end
    %                 end
    %             else
    %                 CursorO(:,tt) = CursorO(:,tt-1);
    %             end
    %             tindex = tindex+8;
    %         end
    %
    %         coordToState = 40.95;  %This is from the DoFeedback func. of DasherSpeller.cpp
    %         CursorO = fix(CursorO.*coordToState);
    %         %Only interested in the second index - the Y position of the cursor
    %         CursorO = (CursorO(1,:));
    %         CursorO = repmat(CursorO,8,1);
    %         Output4 = reshape(CursorO,1,size(CursorO,1)*size(CursorO,2));
    %
    %
    %         plot(Output4,'Color',colors{cv})
    %     end
    
end

% clear OutPutC CumOutPut
% for tt = 1:size(P,3)
%     outstep =0;
%     for cc = 1:size(ClassifierMatrix,1)
%         outstep = outstep + ClassifierMatrix(cc,4)*...
%             P(ClassifierMatrix(cc,1),ClassifierMatrix(cc,2),tt);
%     end
%     OutPutC(tt) = outstep;
% end
%
% % OutPutC = OutPutC + .1;
%
% for tt = 1:size(P,3)
%     if TC(featinds(tt)) == 0
%         CumOutPut(tt) = 0;
%     else
%         CumOutPut(tt) = CumOutPut(tt-1)+OutPutC(tt);
%     end
% end
% OutPutC = repmat(OutPutC,8,1);
% OutPutC = OutPutC(:);
% CumOutPut = repmat(CumOutPut,8,1);
% CumOutPut = CumOutPut(:);
%
% [bc ac] = butter(5,5/(fs/2),'low')
% OutPutC2 = filtfilt(bc,ac,OutPutC);
%
% figure
% plot(OutPutC);
% hold on
% plot(OutPutC2,'g')
% plot(100*TC,'r');
% plot(CumOutPut,'k');

if Figures_On==1
    if L_trialend(end)>length(Output4)
        L_trialend(end) = length(Output4);
    end
    if R_trialend(end)>length(Output4)
        R_trialend(end) = length(Output4);
    end
    
    plotAcc = Output4([L_trialend; R_trialend])'-2047;
    acc_fromPlot(plotAcc<0) = TwoClass(1);
    acc_fromPlot(plotAcc>=0) = TwoClass(2);
    
    acc_fromPlot = sum(acc_fromPlot'==classAvg)/...
        length(classAvg)
end

cursAf(cursA>NormalizerOffset) = TwoClass(2);
cursAf(cursA<=NormalizerOffset) = TwoClass(1);

acc_fromCfier = sum(cursAf==classAvg')/length(classAvg)




sprintf(['Accuracies may not be the same because acc_fromCfier is \n'...
    'calculated from a fixed (3s) window, where the one calculated from\n'...
    'the plot is may be shorter if the target was reached sooner.'])

end
%Files updated June/July 2015 for use with second round of ALS patients
%In addition to some minor changes, added the ability to use the RSVP
%speller, as well as a regularized (LASSO) regression classifier option
%(Cfier == 4).

% for v2, implemented classification with a range of sequences.  Repurposed
% the SequencestoRemove variable into the classification scheme, so that we
% can see how the classification accuracy changes with removal of
% sequences. 9/17/15
%
% in v3, Added the ability to remove epochs from teh data, while retaining
% the ability to simulate using less sequences for classification under
% both the non-averaged and averaged ERP scenarios (Avv=0 and Avv=1).
% 10/2/15
%
% 12/15/15 - Modified v3 to allow for runs from multiple sessions to
% beGoalLll
% used in classifier
%
% ~4/1/16 created BuildClassifier_teleBCI.m
%
% 4/20/16 - Changed code to allow for use of accumulate evidence.  Deals
% with the ability to aggregate flashes from multiple rounds of flashes
% 6/1/16 - teleBCI2 takes the location of the data as an input
%
% 2/3/17 - teleBCI3 removes the restriction that input runs have to have
% the same number of sequences to be processed.  This required changes to
% the way the log file is read, eye tracking data is processed, and data is
% accumulated.
%
% 5/11/17 - teleBCI4 reorganization of the data organization steps.  It was
% slow when combining files with disparate sequence numbers.
%
% 5/24/17 - teleBCI5 changes to audio speller.  When defining target vs
% non-target sequences, i think i should use the attended target as choice
% and the non attended target as the not choice.  Do not include the beeps
% in the classifier at all.
%
% 6/27/17 - Changes to heartbeat code


function [CM, allsums, targetstim, WORDD_full, skiptest, TTSpell, ChosenLetterF,...
    TrainD,ClassLabel, ChanLabel, SC, SR, LDAtimes, prmfileshort, eyye, range, NotDownD] = ...
    BuildClassifier_teleBCI7(Name,Sessions,Runs,BuildClassifier,Cfier,...
    Batch,RerefVal,Art,DoCSP,CSession,CRun, Figures_On, bcifold, ctype)

if Figures_On
    %    close all;
else
end
wtbrr = waitbar(0,'Creating Classifier...');
rng('shuffle')


disp('MAKE SURE ALIGN CHANNELS IS 0')
disp('When recording, data is saved after being filtered, ')
disp('no other processing has gone into saved data');
%Load Data and  Log File

cdp = cd;

if ~exist('bcifold','var')
    [~, sl] = regexp(cdp,'\Users\');
    se = regexp(cdp,'\');
    tmp = find(sl==se);
    pcusr = cdp(sl+1:se(tmp+1)-1);
    bcifold = ['C:\Users\' pcusr '\Documents\BCI2000_5300'];
end


Avv = 0; %Classify from individual trials(0) or averages of stimulus codes(1)
% Figures_On = 0;
HBF = 1;

if Batch == 0
    RerefVal  = 0;  %0 = no reref, 1 = CAR, 2 = Lap, 3 = Bipolar (Neither Lap nor Bip work well)
    Art=1; %Artifact Data (0 - nothing, 1 - remove regions of artifact, 2 - regression)
    DoCSP = 0;
    Figures_On = 1;
    
end

if Cfier == 1
    ReducT = 3; %No data reduction for stepwise regression
else
    ReducT = 3; %Perform data reduction for LDA - 1=decimate, 2=signif test, 3 = Moving average,downsample
end

SequencestoRemove = 0;  %For 20 flashes, there are 10 sequences.


if BuildClassifier == 0
    if Batch == 0;
        CSession = input('What Session? (input as string)   ');
        CRun = input('What Run? (input as string)   ');
    elseif Batch == 1;
        CRun = cell2mat(cat(2,CRun{:}));
        CSession = cell2mat(CSession);
    end
    %     NameBase = ['C:\Documents and Settings\amg5106\My Documents\'...
    %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' Name CSession '\'...
    %         Name 'S' CSession 'R' CRun '_'];
    NameBase = [bcifold '\data\' Name CSession...
        '\' Name 'S' CSession 'R' CRun '_'];
    
    NameBase2 = ['P:\ALS Proj Data\' Name CSession '\'...
        Name 'S' CSession 'R' CRun '_'];
end



%% Load Data
tic
Data = [];
StimulusCode = [];
StimulusType = [];
InputText = [];
SourceTime = [];
StimulusTime = [];
SelectedStimulus = []; % For audio speller
NewRunName = [];
TTSpell = [];
GazeX = [];
GazeY = [];
SequencePhase = [];
AllTrials = [];
AppPos = [];

kk = 1; max_sess = 0;
for ss = 1:length(Sessions)
    Session = Sessions{ss};
    curr_sess = str2double(Session);
    if curr_sess>max_sess
        max_sess = curr_sess;
        max_sess_i = ss;
    end
    for rr = 1:length(Runs{ss})
        %     [a b c d] = load_bcidat(['C:\Documents and Settings\amg5106\My Documents\'...
        %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' ...
        %         Name Session '\' Name 'S' Session 'R' Run{rr} '.dat'],'-calibrated');
        try %try first to find the data on the computer
            [a b c d] = load_bcidat([bcifold '\data\' ...
                Name Session '\' Name 'S' Session 'R' Runs{ss}{rr} '.dat'],'-calibrated');
            disp(['Time to load data ' num2str(toc)]);
        catch %If that doesnt work, try on the cluster server space
            [a b c d] = load_bcidat(['P:\ALS Proj Data\' ...
                Name Session '\' Name 'S' Session 'R' Runs{ss}{rr} '.dat'],'-calibrated');
            disp(['Time to load data ' num2str(toc)]);
        end
        
        
        if isfield(c,'SineChannelX') %overwrite headset type if Signal Generator
            capt = '_siggen';
        end
        
        
        if size(a,2)  == 14 %Emotiv Headset
            capt = 'emotiv';
            ChannelNames = cell(14,1);
            EEGloc = 1:14; EOGloc = [];
        elseif size(a,2)  == 22 %standard EEG cap (2 gusb)
            capt = '22_3EOG_2amp';
            ChannelNames = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8',...
                'P7','P3','Pz','P4','P8','O1','O2'};
            EEGloc = 1:19; EOGloc = 20:22;
        elseif size(a,2)  == 16 % 16 chan gNautilus
            capt = 'gNautilus16';
            ChannelNames = {'Fp1','Fp2','F3','Fz','F4','T7','C3','Cz','C4','T8',...
                'P3','Pz','P4','PO7','PO8','Oz'};
            EEGloc = 1:16; EOGloc = 1;
        elseif size(a,2)  == 8 % 8 chan gNautilus or Signal Generator
            if isfield(c,'DeviceIDs')
                capt = 'gnautilus_8wet';
                EEGloc = 1:8;
                if sum(strcmp(Name,{'T02','T05'}))
                    ChannelNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};
                    EOGloc = 1;
                else
                    disp('DEFAULTING TO NEWER FLEXCAP CONFIG - FZ = CH2')
                    ChannelNames = {'Cz','Fz','P3','Pz','P4','PO7','PO8','Oz'};
                    EOGloc =2;
                end
            else
                capt = 'siggen';
                ChannelNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};
                EEGloc = 1:8; EOGloc = 1;
            end
            
        elseif size(a,2)  == 3
            capt = '3ch';
            ChannelNames = cell(3,1);
        else
            disp('Undefined electrode configuration!')
            return
        end
        
        
        if strcmp(capt,'siggen')
            Art = 0;
        end
        
        
        
        
        
        
        
        
        % Shift Data? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %     kk = 1;
        %     for jj = [16 17]
        %         [MM, ~] = peakdet(a(:,jj),50);
        %         MaxTab{kk} = MM(:,1);
        %         kk = kk+1;
        %     end
        %
        %     %Positive differences are when amp2 leads amp1.
        %     CumDiff = []; CumLoc = [];
        %     for jj = 1:length(MaxTab{1})
        %         CumDiff = [CumDiff; MaxTab{1}(jj)-MaxTab{2}];
        %         CumLoc = [CumLoc; repmat(MaxTab{1}(jj),length(MaxTab{2}),1)];
        %     end
        %
        %     keepem = CumDiff<10&CumDiff>-10;
        %     CumDiff = CumDiff(keepem);
        %     CumLoc = CumLoc(keepem);
        %
        %
        %     if Figures_On == 1
        %         figure
        %         tx(1) = subplot(221); plot(a(:,1)); hold on; plot(a(:,21),'r');
        %         subplot(2,2,[2 4]); hist(CumDiff,[-10:10]); xlim([-10 10])
        %     end
        %
        %     Mshift = mode(CumDiff);
        %     if Mshift > -2 && Mshift < 2
        %     a(:,17:end) = circshift(a(:,17:end),[Mshift 0]);
        %     disp(['SHIFTED by' num2str(Mshift)])
        %     end
        %
        %     if Figures_On == 1
        %         tx(2) = subplot(223); plot(a(:,1)); hold on; plot(a(:,21),'r');
        %         title(['Shift ' num2str(Mshift)]);
        %         linkaxes(tx,'x')
        %     end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp([num2str(c.PreRunDuration.NumericValue) ' ' num2str(c.PostRunDuration.NumericValue) ...
            ' ' num2str(c.PreSequenceDuration.NumericValue) ' ' num2str(c.PostSequenceDuration.NumericValue)]);
        
        fs = c.SamplingRate.NumericValue;
        %Trim Data and StimulusCodes
        trimMax = round((c.PreRunDuration.NumericValue+c.PreSequenceDuration.NumericValue-.5)*fs);
        Data = [Data;  a(trimMax:end,:)];
        StimulusCode = [StimulusCode; b.StimulusCode(trimMax:end)];
        StimulusType = [StimulusType; b.StimulusType(trimMax:end)];
        SourceTime = [SourceTime; b.SourceTime(trimMax:end)];
        StimulusTime = [StimulusTime; b.StimulusTime(trimMax:end)];
        if isfield(c,'ToBeCopied') %audio speller
            TTSpell = [TTSpell c.ToBeCopied.NumericValue'];
            TTSpell_cell{kk} = c.ToBeCopied.NumericValue;
            SelectedStimulus = [SelectedStimulus; b.SelectedStimulus(trimMax:end)];
        else
            TTSpell = [TTSpell cell2mat(c.TextToSpell.Value)];
            TTSpell_cell{kk} = cell2mat(c.TextToSpell.Value);
        end
        if isfield(b,'EyetrackerLeftEyeGazeX')
            GazeX = [GazeX; b.EyetrackerLeftEyeGazeX(trimMax:end)];
            GazeY = [GazeY; b.EyetrackerLeftEyeGazeY(trimMax:end)];
        else
            GazeX = [GazeX; NaN(size(a(trimMax:end,:),1),1)];
            GazeY = [GazeY; NaN(size(a(trimMax:end,:),1),1)];
        end
        SequencePhase = [SequencePhase; b.PhaseInSequence(trimMax:end)];
        NewRunName = strcat(NewRunName,Runs{ss}{rr})
        %         if max(StimulusCode)==16  %Is this the multiflash grid (this) or single flash (next)?
        %             %                 fid = fopen(['C:\Documents and Settings\amg5106\My Documents\'...
        %             %                     'Year2\Research\BCI2000\BCI2000sourcetree\data\' Name Session '\'...
        %             %                     Name 'S' Session 'R' num2str(str2num(Run{rr})) '+_mylogfile.txt'],'r');
        %
        %
        %             fid = fopen([bcifold '\data\' Name Session '\'...
        %                 Name 'S' Session 'R' num2str(str2num(Runs{ss}{rr})) '+_mylogfile.txt'],'r');
        %             if fid==-1
        %                 fid = fopen(['P:\ALS Proj Data\' Name Session '\'...
        %                     Name 'S' Session 'R' num2str(str2num(Runs{ss}{rr})) '+_mylogfile.txt'],'r');
        %             end
        %
        %             itt=textscan(fid,'%s','delimiter','\n');
        %             InputText = [InputText; itt{1}];
        %             Speller = 1; %multiflash grid
        %             fclose(fid);
        %         elseif max(StimulusCode)==12 %The unmodified (RC) P300 speller
        %             Speller = 3;
        %         end
        
        if isfield(c,'ToBeCopied') %audio speller
            Speller = 4;
        else
            %Checkerboard Speller
            if max(StimulusCode) > (c.NumMatrixRows.NumericValue+c.NumMatrixColumns.NumericValue)
                fid = fopen([bcifold '\data\' Name Session '\'...
                    Name 'S' Session 'R' num2str(str2num(Runs{ss}{rr})) '+_mylogfile.txt'],'r');
                if fid==-1
                    fid = fopen(['P:\ALS Proj Data\' Name  '\'...
                        Name 'S' Session 'R' num2str(str2num(Runs{ss}{rr})) '+_mylogfile.txt'],'r');
                end
                
                itt=textscan(fid,'%s','delimiter','\n');
                InputText = [InputText; itt{1}];
                Speller = 1; %multiflash grid
                fclose(fid);
            end
            %RC Speller
            if max(StimulusCode) == (c.NumMatrixRows.NumericValue+c.NumMatrixColumns.NumericValue)
                Speller = 3;
            end
            %Also check for an RSVP speller, which can have any number of stimulus
            %codes, but will always have NumMatrixRow = 1. max(StimulusCode) will
            %equal 2*NumColumns
            if c.NumMatrixRows.NumericValue == 1 %This is the RSVP speller
                Speller = 2;
            end
        end
        
        if ~exist('Speller','var')
            disp('Speller format not detected'); return;
        end
        
        
        
        %Associate Target Codes with Targets
        if Speller == 4
            tdtemp = c.Stimuli.Value(1,:);
            TargSymb(kk,:) = tdtemp
        else
            tdtemp = c.TargetDefinitions.Value;
            if size(tdtemp,2)==1 %Multi Menu
                TargSymb(kk,:) = tdtemp{1}(:,2)'; %-- NOTT CORRECT YET
            else
                TargSymb(kk,:) = tdtemp(:,2)';
            end
            SymbSymb = {'Angry%20','Bag%20','Bed%20','Bored%20','Bathroom%20','Carer%20','Cdrink%20',...
                'Clothes%20','Cold%20','Doctor%20','Family%20','Food%20','Glasses%20','Hdrink%20',...
                'Hear%20','Help%20','Hot%20','Spouse%20','Idk%20','Light%20','Medicine%20','Yes%20',...
                'No%20','Nurse%20','Pain%20','Pill%20','Slippers%20','Telephone%20','Tired%20',...
                'Toilet%20','Tv%20','Wheelchair%20','Walker%20','Worried%20','Watch '};
        end
        
        %Number of channels
        NumChans{kk} = size(a,2);
        
        %Number of sequences per trial.  Two times this number is the number of
        %times each stimulus is flashed each trial.
        NumSequences{kk} = c.NumberOfSequences.NumericValue;
        
        %Number of StimulusCodes
        NumStimCodes = max(StimulusCode);
        
        %Duration of the Stimulus State -- This is how long the stimulus is on
        StateDuration = c.StimulusDuration.NumericValue*fs/1000;
        
        %Number of trials
        NumTrials{kk} = fix((sum(b.StimulusCode>=1)/StateDuration)/(NumSequences{kk}*NumStimCodes));
        
        if isfield(b,'EyetrackerLeftEyeGazeX')
            AppPos = [AppPos; repmat([c.WindowLeft.NumericValue c.WindowTop.NumericValue ...
                c.WindowWidth.NumericValue c.WindowHeight.NumericValue],NumTrials{kk},1)];
        end
        
        %Assign each trial to a run number
        AllTrials = [AllTrials; kk*ones(NumTrials{kk},1)];
        
        %Load Spatial Filter
        SpatFiltUsed{kk} = c.SpatialFilter.NumericValue;
        
        %Load Classifier
        ClassifierUsed{kk} = c.Classifier.NumericValue;
        ClassifierUsed{kk}(:,4) = cellfun(@str2num,c.Classifier.Value(:,4));
        
        kk = kk+1;
    end
end
Data = Data';
% disp('DONT FORGET TO REMOVE THISSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS');
% Data = Data*.1;
if ~isempty(strfind(capt,'gnautilus'))
    range = (round(.45*fs):round(1.2*fs)); %Use this data range for classifiation
else
    range = 1:round(.8*fs);
end
disp(['Time to load data ' num2str(toc)])
waitbar(.05,wtbrr);
%  range = (1:round(fs));

%% Check Variables
%Check if free spelling
FreeSpell = sum(StimulusType)==0&sum(StimulusCode)~=0;

%Are channels consistant across runs?
if sum(diff(cell2mat(NumChans))) ~= 0
    disp('Different number of channels in each run!');
end
NumChans = NumChans{1};

% %Do the channel labels match the number of channels in Data?
% if length(EOGloc)+length(EEGloc)~=size(Data,1)
%     disp('Channel labels incorrect'); return
% end

%Are the spatial filters the same same size?
try %Are classifiers the same size
    s2mC = cell2mat(SpatFiltUsed);
    
    if sum(sum(diff(reshape(s2mC,size(s2mC,1),size(s2mC,2)/rr,rr),[],3))) ~= 0
        disp('Spatial Filters have different values in each run!');
    end
    SpatFiltUsed = SpatFiltUsed{1};
catch err
    disp('Spatial Filters  are different sizes in each run!');
end
%If the same size, do they contain the same values?


try %Are classifiers the same size
    c2mC = cell2mat(ClassifierUsed);
    %If the same size, do they contain the same values?
    if sum(sum(diff(reshape(c2mC,size(c2mC,1),size(c2mC,2)/rr,rr),[],3))) ~= 0
        disp('Classifiers have different values in each run!');
    end
    ClassifierUsed = ClassifierUsed{1};
catch err
    disp('Classifiers are different sizes in each run!');
end


% AG 2/3/17 attempting to remove this restriction
% %Are number of sequences the same
% if sum(cell2mat(NumSequences))/(kk-1) ~= NumSequences{1}
%     disp('Not the same number of sequences in each run, code will not run correctly');
% end
% NumSequences = NumSequences{1};
% if SequencestoRemove >= NumSequences
%     disp('Removed too many sequences');
% end

%Mark the start times of each trial
TrialStartTime = find(SequencePhase==1 & [diff(SequencePhase); 0]);

if Speller==4
    NR = 1; NC = 1;
else
    %Number of rows, columns, stimuli
    NR = c.NumMatrixRows.NumericValue;
    NC = c.NumMatrixColumns.NumericValue;
end

% if size(tdtemp,2)==1 %Multi Menu
%     NS = size(tdtemp{1},1);
% else
%     NS = size(tdtemp,1);
% end
%Check state duration again
for scc = 1:NumStimCodes
    sdd = find(StimulusCode == scc);
    kk(scc)=2;
    try
        while((sdd(kk(scc))-sdd(kk(scc)-1))==1)
            kk(scc) = kk(scc)+1;
        end
    end
end
if mean(kk) == kk(1)
    StateDuration2 = mean(kk)-1;
else
    disp('We have a problem');
end

if StateDuration ~= StateDuration2
    disp('Specified stimulus time not consistant with actual on time');
    StateDuration = StateDuration2;
end


%% Check Timing of Data
%Check for Dropped Samples (in the case of BCI2000, the data is repeated in
%the next block if the processing exceeds the roundtrip time)
Timerr = diff(SourceTime(1:c.SampleBlockSize.NumericValue:end));
StTimerr = diff(StimulusTime(1:c.SampleBlockSize.NumericValue:end));
thTime = c.SampleBlockSize.NumericValue/fs*1000;

BadTime = Timerr<thTime-5 | Timerr>thTime+5;
if Figures_On == 1
    figure('Name','Stimulus time and Source time and Stimulus Code')
    plot(StTimerr); hold on;
    plot(Timerr,'r'); ylim([-12 250]);
    plot(StimulusCode(1:c.SampleBlockSize.NumericValue:end),'g');
    plot(-BadTime*10,'k')
end
%figure
%plot(Data(:,1)*1000); hold on; plot(b.SourceTime,'r');



%% Artifact Correction
StimulusCodeUNArt = StimulusCode;
StimulusTypeUNArt = StimulusType;



%Remove dead data
SmInt = [];
for i= 256:256:size(Data,2)-255
    if std(Data(EOGloc(1),i:i+255))<1;
        SmInt  = [SmInt i];
    end
end
if ~isempty(SmInt)
    Srange = fix(-3*fs:5*fs);
    SmInt = repmat(SmInt,length(Srange),1)+repmat(Srange,length(SmInt),1)';
    SmInt = SmInt(:);
    SmInt = SmInt(SmInt>0 & SmInt<size(Data,2));
    SmLoc = zeros(size(Data,2),1);
    SmLoc(SmInt)=1;
    StimulusType(logical(SmLoc))=0; %= StimulusType(~ArtLoc);
    StimulusCode(logical(SmLoc))=0; %= StimulusCode(~ArtLoc);
else
    SmLoc = [];
end

%Remove superlarge data artifacts 1s before, 5 seconds after
LInt = find(abs(Data(EOGloc(1),:))>1e3);
if ~isempty(LInt)
    Lrange = fix(-1*fs:5*fs);
    LInt = repmat(LInt,length(Lrange),1)+repmat(Lrange,length(LInt),1)';
    LInt = LInt(:);
    LInt = LInt(LInt>0 & LInt<size(Data,2));
    LLoc = zeros(size(Data,2),1);
    LLoc(LInt)=1;
    StimulusType(logical(LLoc))=0; %= StimulusType(~ArtLoc);
    StimulusCode(logical(LLoc))=0; %= StimulusCode(~ArtLoc);
else
    LLoc = [];
end

%Ocular artifacts
if isempty(EOGloc)
    disp('Cannot reject artifacts because no EOG channels were recorded');
    Art = 0;
end
if Art == 1 %Remove data displaying artifacts
    disp('Removing ocular artifacts....');
    
    ArtInt = find(Data(EOGloc(1),:)>50&circshift(Data(EOGloc(1),:),[0 1])<50);
    ArtInt = [ArtInt find(Data(EOGloc(1),:)<-50&circshift(Data(EOGloc(1),:),[0 1])>-50)];
    
    Arange = fix(-1*fs:.5*fs);  %This was changed 5/29/17 from -.1 to 1
    %because (a) we have to account for the .45 second lag in EEG data
    %to account for the fact that the large peak was corrupting the EEG
    %data for preceding stimulus codes that were only set to zero .1
    %seconds before.
    ArtInt = repmat(ArtInt,length(Arange),1)+repmat(Arange,length(ArtInt),1)';
    %     ArtInt = ArtInt';
    ArtInt = ArtInt(:);
    ArtInt = ArtInt(ArtInt>0 & ArtInt<size(Data,2));
    ArtLoc = zeros(size(Data,2),1);
    ArtLoc(ArtInt)=1;
    
    StimulusType(logical(ArtLoc))=0; %= StimulusType(~ArtLoc);
    StimulusCode(logical(ArtLoc))=0; %= StimulusCode(~ArtLoc);
    
    if Figures_On == 1
        figure('Name','EOG channel 1, locations of artifact, and stimulus code')
        plot(Data(EOGloc(1),:)); hold on;
        plot(50*ArtLoc,'r');
        plot(50*SmLoc,'c');
        plot(50*LLoc,'m');
        %     figure
        %     tt(1) = subplot(211); plot(Data(EEGloc,:)');
        %     tt(2) = subplot(212); plot(Data(EOGloc,:)');
        %     linkaxes(tt)
    end
    
    %Find incomplete trials
    sdtemp = find(diff(double(StimulusCode))~=0);
    nftemp = diff(sdtemp)<StateDuration;
    nfstart = sdtemp(nftemp)+1;
    nfend = sdtemp(find(nftemp)+1);
    
    for ff = 1:length(nfstart)
        StimulusCode(nfstart(ff):nfend(ff)) = 0;
        StimulusType(nfstart(ff):nfend(ff)) = 0;
    end
    
    
    %     ismember(diff(double(StimulusCode))>0
    
    
    if Figures_On == 1
        plot(StimulusCode,'k');
    end
    
    ArtWeights = eye(size(Data,1));
    ArtWeights2 = eye(length(EEGloc),size(Data,1));
    NumChans = length(EEGloc);
    disp(['Done with Art ' num2str(toc)]);
elseif Art == 2 %Artifact regression
    if BuildClassifier == 1
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
            if length(EOGloc)==3
                ArtInt = find(Data(EOGloc(1),:)-Data(EOGloc(2),:)>75 | ...
                    Data(EOGloc(1),:)-Data(EOGloc(2),:)<-75 | ...
                    Data(EOGloc(1),:)-Data(EOGloc(3),:)>75 | ...
                    Data(EOGloc(1),:)-Data(EOGloc(3),:)<-75);
            else
                ArtInt = find(Data(EOGloc(1),:)>75 | ...
                    Data(EOGloc(1),:)<-75);
            end
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
            if length(EOGloc)==3
                figure('Name','Differential EOG and location of artifacts')
                plot(Data(EOGloc(1),:)-Data(EOGloc(2),:))
                hold on
                plot(Data(EOGloc(1),:)-Data(EOGloc(3),:))
                plot(100*ArtLoc,'r')
            else
                figure('Name','EOG and location of artifacts')
                plot(Data(EOGloc(1),:)); hold on;
                plot(100*ArtLoc,'r');
            end
        end
        
        EOGD = Data(:,ArtLoc);
        
        if ~isempty(EOGD)
            %the function covm adds an additional column of ones in front of the data
            %and is necessary for regress_eog.m
            if length(EOGloc)==3
                [R] = regress_eog(covm(EOGD','E'),EEGloc, ...
                    sparse([EOGloc(1),EOGloc(3),EOGloc(2),EOGloc(1)],[1,1,2,2],[1,-1,1,-1]));
            elseif length(EOGloc)==2
                [R] = regress_eog(covm(EOGD','E'),EEGloc, ...
                    sparse([EOGloc(1),EOGloc(2)],[1,1],[1,-1]));
            else
                [R] = regress_eog(covm(EOGD','E'),EEGloc, EOGloc);
            end
            %Create full matrix for online artifact reduction
            %I believe this is the way they say to do it (pad Data with a channel
            %of ones -- this introduces a bias to the output channel) (see DD2 below).
            %However, this padding is not something I want to do online, and since
            %it is only a bias, we can remove the first column of ArtWeights.
            ArtWeights = full(R.r0)';
            ArtWeights2 = ArtWeights(EEGloc,:);
            %         ArtWeights = full(R.r1)';
            %         ArtWeights2 = ArtWeights(:,2:end);
            %DD2 = [ones(size(Data,1),1),Data] * ArtWeights';
        else
            ArtWeights2 = eye(length(EEGloc),size(Data,1));
        end
        
        NumChans = length(EEGloc);
        
        
    elseif BuildClassifier == 0
    end
    disp(['Done with Art ' num2str(toc)]);
else
    ArtWeights = eye(size(Data,1));
    ArtWeights2 = eye(length(EEGloc),size(Data,1));
    NumChans = length(EEGloc);
    %     StimulusCodeUNArt = StimulusCode;
    %     StimulusTypeUNArt = StimulusType;
    disp('No artifacting done');
end
%Using the correction coefficients, transform entire training run to reduce
%artifact
%
if Figures_On == 1
    if BuildClassifier == 1
        Dw = ArtWeights2*Data;
        %
        %
        %Plot the F and O raw and artifacted, as well as EOG channels
        figure('Name','Original and Corrected Ch1 data, Ch2 data, and EOG channels')
        tt(1)=subplot(311); plot(Data(1,:)); hold on; plot(Dw(1,:),'r');
        tt(2)=subplot(312); plot(Data(2,:)); hold on;
        plot(Dw(2,:),'r');
        tt(3)=subplot(313); plot(Data(EOGloc,:)'); hold on
        linkaxes(tt,'x');
        
        clear Dw;
    end
end
%
% %Plot the spectra of the raw vs rerefed and artifacted channels
%     figure
%     kk=1;
%     for ich =[EEGloc(1) EEGloc(2) EEGloc(end-1-length(EOGloc)) EEGloc(end-length(EOGloc))]
%         subplot(2,2,kk);
%         [Sp1,ff1]=pwelch(Data(:,ich),2*fs,fs,2*fs,fs); hold on;
%         plot(ff1,10*log10(Sp1));
%         [Sp2,ff2]=pwelch(DD(:,ich),2*fs,fs,2*fs,fs);
%         plot(ff2,10*log10(Sp2),'r');
%         kk=kk+1;
%     end
%
% DD = DD(:,EEGloc);
% NumChans = length(EEGloc);






%% Heartbeat filter

%Typically the heartbeat shows up in PO7, PO8, and Oz, but it can
%corrupt all channels.  Find the heartbeat in these three channels,
%by filtering out low frequency data, performing peak detection, creating a
%matched filter, refining the heartbeat selections, and subtracting the
%heartbeat template.
if HBF == 1
    [bb,aa] = butter(4,[4 60]/(fs/2));    
    HBdata = filtfilt(bb,aa,Data')';
    temp = mean(HBdata(6:8,:),1);
    
    
    %Expect roughly one beat per second
    expectedbeats = length(temp)/(fs);
    i = 8;
    maxt = [];
    while size(maxt,1)<expectedbeats
        if i == 3
            disp('no heartbeat detected')
            break
        end
        [maxt, mint] = peakdet(temp,i*std(temp));
        i=i-1;
    end
    
    
    
    if ~isempty(maxt)
        maxt(maxt(:,1)<fs|maxt(:,1)>size(Data,2)-fs,:) = [];
        bpm = fs*60./diff(maxt(:,1));
        %remove heartbeats greater than 200 bpm
        bpmi = bpm>200 | circshift(bpm>200,1);
        keepbeep = true(length(bpmi),1);
        counter = []; start = 0;
        for jj = 1:length(bpmi)
            if bpmi(jj) == 0
                if start == 1
                [xi,yi] = sort(counter(:,2));
                keepbeep(counter(yi(1:end-1),1)) = 0;
                counter = []; start = 0;
                end
            else
                counter = [counter; [jj maxt(jj,2)]];
                start = 1;
            end
        end
        if ~isempty(counter)
               [xi,yi] = sort(counter(:,2));
                keepbeep(counter(yi(1:end-1),1)) = 0;
                counter = []; start = 0;
        end
        maxt = maxt(find(keepbeep),:)
        bpm2 = fs*60./diff(maxt(:,1));
        statschan = [mean(bpm2) std(bpm2)];
        maxval = maxt(:,2);
        maxloc = maxt(:,1);
        if Figures_On==1
            figure
            yy(1) = subplot(3,4,1:3); plot(temp/max(temp),'Color',[.6 .6 .6]); hold on;
            scatter(maxloc, maxval/max(temp),100,'.k');
            subplot(3,4,4); hist(bpm2,200);
        end
    else
        statschan = [0 100];
    end
    
    if statschan(2)<40
        disp(['Preliminary heartbeat at ' num2str(statschan(1)) ' bpm, creating matched filter'])
        keepbeep2 = bpm2>statschan(2)-20 & bpm2<statschan(1)+20;
        tempind = repmat(maxloc(keepbeep2,1),1,fs)+repmat(-(fs/2-1):fs/2,sum(keepbeep2),1);
        template = nanmean(temp(tempind));
        
        %         temp2 = conv(temp,template,'same');
        temp3 = conv(temp,fliplr(template),'same');
        
        %         figure
        %         xx(1) = subplot(311); plot(temp);
        %         xx(2) = subplot(312); plot(temp2);
        %         xx(3) = subplot(313); plot(temp3);
        %         linkaxes(xx,'x');
        
        i = 6;
        maxt3 = [];
        while size(maxt3,1)<expectedbeats
            [maxt3, mint3] = peakdet(temp3,i*std(temp3));
            i=i-1;
            if i == 0
                
                break
            end
        end
        if ~isempty(maxt3)
            maxt3(maxt3(:,1)<fs|maxt3(:,1)>size(Data,2)-fs,:) = [];
            bpm3 = fs*60./diff(maxt3(:,1));
            %remove heartbeats greater than 200 bpm
            bpm3i = bpm3>200 | circshift(bpm3>200,1)
            keepbeep3 = true(length(bpm3i),1);
            counter = []; start = 0;
            for jj = 1:length(bpm3i)
                if bpm3i(jj) == 0
                    if start == 1
                        [xi,yi] = sort(counter(:,2));
                        keepbeep3(counter(yi(1:end-1),1)) = 0;
                        counter = []; start = 0;
                    end
                else
                    counter = [counter; [jj maxt3(jj,2)]];
                    start = 1;
                end
            end
            if ~isempty(counter)
                [xi,yi] = sort(counter(:,2));
                keepbeep3(counter(yi(1:end-1),1)) = 0;
                counter = []; start = 0;
            end
            maxt3 = maxt3(find(keepbeep3),:)
            bpm4 = fs*60./diff(maxt3(:,1));
            statschan3 = [mean(bpm4) std(bpm4)];
            maxval3 = maxt3(:,2);
            maxloc3 = maxt3(:,1);
            if Figures_On==1
                yy(2) = subplot(3,4,5:7); hold on; plot(temp3/max(temp3),'Color',[1 .5 .5]);
                scatter(maxloc3, temp3(maxloc3)/max(temp3),100,'.r');
                subplot(3,4,8); plot(template,'k');
            end
            
            Data2 = Data;
            tempind = repmat(maxt3(:,1),1,fs)+repmat(ceil(-fs/2+1:fs/2),size(maxt3,1),1);
            wnd = window(@hanning,fs);
            for i = 1:8
                tempd = reshape(HBdata(i,tempind'),fs,size(tempind,1));
                %remove artifactual epochs
                bade = var(tempd)>5*abs(mean(var(tempd)));
                tempd(:,bade) = [];
                signature(:,i) = wnd.*mean(tempd,2);
                
                for j = 1:size(maxt3,1)
                    Data2(i,maxt3(j,1)+ceil(-fs/2+1:fs/2)) = Data2(i,maxt3(j,1)+ceil(-fs/2+1:fs/2))-signature(:,i)';
                end
            end
            
            if Figures_On==1
                subplot(3,4,12); plot(signature);
                yy(3) = subplot(3,4,9:11);  plot(mean(Data(6:8,:),1),'k'); hold on; 
                plot(mean(Data2(6:8,:),1),'Color',[.5 .5 .5]);
                plot(mean(Data(1,:),1)-50,'b'); hold on; plot(mean(Data2(1,:),1)-50,'Color',[.5 .5 1]);
                linkaxes(yy,'x');
            end
            
        else
            disp('no heartbeat filter created');
        end
        
        
        
        %         %I dont think ICA is ideal for this situation.  First, it takes
        %         %long.  Second there are only eight channels.
        %         [weights, sphere] = runica(Data);
        %         temp = weights*Data;
        %         figure
        %         plot(temp');
        %
        %         for i = 1:8
        %             [maxt, mint] = peakdet(temp(i,:),5*std(temp(i,:)));
        %             bpm = fs*60./diff(maxt(:,1));
        %             %remove heartbeats greater than 100 bpm
        %             keepbeep = bpm<=100;
        %             bpm(~keepbeep,:) = [];
        %             scatter(maxt(keepbeep,1),maxt(keepbeep,2));
        %             subplot(1,4,4); hist(bpm,200);
        %             statsica(i,:) = [mean(bpm) std(bpm)];
        %         end
        %         [~, candidate] = min(statsica(:,2));
        %         if abs((statsica(candidate,1)-statschan(1))/statschan(1)*100)<5 %bpm
        %             %from the candidate ica channel is less than 5% off the bpm in the
        %             %PO/O channels.  Remove this signal from data.
        %             invweights = inv(weights);
        %             invweights(:,candidate) = 0;
        %             temp2 = invweights*temp;
        
        
        if Figures_On == 1
            figure
            xx(1) = subplot(211); plot(Data'); hold on; plot(StimulusCodeUNArt,'k','LineWidth',1.5);
            xx(2) = subplot(212); plot(Data2'); hold on; plot(StimulusCode,'k','LineWidth',1.5);
            linkaxes(xx);
        end
        Data = Data2;
    else
        disp('Heartbeat signal too noisy... skipping filter');
    end
else
    disp('Heartbeat filter disabled')
end


clear HBdata



%% Rereference Data
switch RerefVal
    case 0
        RefFilt = eye(length(EEGloc));
        %         DataF = Data*SpatFilt;
    case 1
        RefFilt = (-1/length(EEGloc))*(ones(length(EEGloc))-eye(length(EEGloc)));
        RefFilt = RefFilt + eye(length(EEGloc));
        %         SpatFilt = blkdiag(SpatFilt, eye(length(EOGloc)));
        %         DataF = Data*SpatFilt;
    case 2
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
        NumChans_old = NumChans;
        NumChans = size(RefFilt,1);
    case 3
        if size(Data,1)==16
            RefFilt = [0 0 0 0 -1 0 0 0 0 0 0 0 1 0;
                0 0 0 0 -1 0 0 0 0 0 0 0 0 1];
        elseif size(Data,1)==22
            RefFilt = [0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 1 0;
                0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 1];
        end
        NumChans_old = NumChans;
        NumChans = 2;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find code stimuli associations (contained in SC and SR), as well as the
%stimuli sequences for each trial.



if Speller == 1  %Use the myfile.txt associated with CB speller runs to
    %find stimulus codes and target associations
    
    for i = 1:length(InputText)
        RunEnd(i) = ~isempty(strfind(InputText{i},'******************'));
    end
    Runind = find(RunEnd==1);
    SCindF = []; SRindF = []; NEXTmindF = []; GoalLindF = [];
    ChosenLetterF = []; GoalLetterF = []; INITmindF = []; GoalTrialF = 0;
    NEXTgoalletterF = [];
    for rr = 1:length(Runind);
        if rr == 1;
            RunRange = 1:Runind(1);
        else
            RunRange = Runind(rr-1):Runind(rr);
        end
        clear Begin SCind SRind INITmind NEXTmind ChosenLind GoalLind
        for ll = RunRange
            GoalLind(ll) = ~isempty(strfind(InputText{ll},'Goal Text,'));
            ChosenLind(ll) = ~isempty(strfind(InputText{ll},'Selected Text'));
            Begin(ll) = ~isempty(strfind(InputText{ll},'*****START'));
            SCind(ll) = ~isempty(strfind(InputText{ll},'SC'));
            SRind(ll) = ~isempty(strfind(InputText{ll},'SR'));
            INITmind(ll) = ~isempty(strfind(InputText{ll},'INIT_mSequence'));
            NEXTmind(ll) = ~isempty(strfind(InputText{ll},'NEXT_mSequence'));
        end
        %Ignore the last SC, SR, NEXTmind, and 2 GoalLetters
        GoalLind2 = find(GoalLind);
        ChosenLind2 = find(ChosenLind)+1; %this is the distance between the selected letter and the repeat of the next goal letter.
        trialkeep = ~ismember(GoalLind2,ChosenLind2);
        trialkeep(end)=0;
        GoalLindF = [GoalLindF (GoalLind2(trialkeep))];
        GoalLetter = InputText(GoalLind2(trialkeep));
        
        [GoalTrialbeg GoalTrialend] = cellfun(@(x) regexp(x,'\s\d*\s='),GoalLetter,'UniformOutput',false);
        GoalTrial = cellfun(@(x,y,z) str2num(x(y+1:z-2)),GoalLetter,GoalTrialbeg,GoalTrialend);
        stind = regexp(GoalLetter,'('); endind = regexp(GoalLetter,')');
        for xy = 1:length(stind)
            GoalLetter{xy} = GoalLetter{xy}(stind{xy}+1:endind{xy}-1);
        end
        
        
        GoalTrial = GoalTrial+max(GoalTrialF);
        GoalTrialF = [GoalTrialF; GoalTrial];%This
        %vector is used to reconcile the number of trials and the number of letters spelled
        %             GoalLetter = TTSpell_cell{rr}(GoalTrial)';
        GoalLetterF = [GoalLetterF; GoalLetter];
        ChosenLetter = InputText(ChosenLind(1:end-1));
        ChosenLetter = cellfun(@(x) x(18:end),ChosenLetter,'UniformOutput',false);
        
        %             ChosenLetter = ChosenLetter(:,18);
        ChosenLetterF = [ChosenLetterF; ChosenLetter];
        Begind = find(Begin==1);
        SCind = find(SCind==1)';
        FSCi = find(SCind < Begind);
        SCindF = [SCindF; SCind(FSCi(end):end-1)];
        SRind = find(SRind==1)';
        FSRi = find(SRind < Begind);
        SRindF = [SRindF; SRind(FSRi(end):end-1)];
        INITmind = find(INITmind==1)';
        INITmindF = [INITmindF; INITmind(1:end-1)];
        NEXTmind = find(NEXTmind==1)';
        omitN = NumSequences{rr}:NumSequences{rr}:NumTrials{rr}*NumSequences{rr};
        NEXTmind(omitN) = [];
        NEXTmindF = [NEXTmindF; NEXTmind];
        tmpind = repmat(GoalTrial,1,NumSequences{rr})';
        tmpind(omitN) = [];
        NEXTgoalletterF = [NEXTgoalletterF; tmpind(:)];
    end
    GoalTrialF = GoalTrialF(2:end);
    
    
    
    
    %         for j = 1:NumSequences-1
    %             tempNEXTm = InputText{NEXTmindF((NumSequences-1)*(i-1)+j)}(16:end);
    %             tempNEXTm = regexp(tempNEXTm,' ','split')'; tempNEXTm(cellfun(@isempty,tempNEXTm)) = [];
    %             NEXTm{j,i} = tempNEXTm;
    %         end
    %     end
    %     if isempty(NEXTm)
    %         mSeq = INITm;
    %     else
    %         mSeq = cat(1,INITm,NEXTm);
    %     end
    
    
    
    
    NEXTm = cell(1,length(GoalLetterF));
    for i = 1:length(GoalLetterF)
        tempSC = InputText{SCindF(i)}(4:end);
        tempSC = regexp(tempSC,' ','split')'; tempSC(cellfun(@isempty,tempSC)) = [];
        SC{i} = tempSC;
        tempSR = InputText{SRindF(i)}(4:end);
        tempSR = regexp(tempSR,' ','split')'; tempSR(cellfun(@isempty,tempSR)) = [];
        SR{i} = tempSR;
        tempINITm = InputText{INITmindF(i)}(16:end);
        tempINITm = regexp(tempINITm,' ','split')';
        tempINITm = cellfun(@(x) str2num(x),(tempINITm(1:end-1,:)));
        INITm{i} = tempINITm;
        
        %tempNEXTm = InputText{NEXTmindF((NumSequences{rr}-1)*(i-1)+j)}(16:end);
        tempNEXTm = InputText(NEXTmindF(NEXTgoalletterF==i));
        tempNEXTm = cellfun(@(x) x(16:end), tempNEXTm,'UniformOutput',false);
        tempNEXTm = regexp(tempNEXTm,'\s','split')';
        tempNEXTm = cat(1,tempNEXTm{:})';
        if ~isempty(tempNEXTm)
            tempNEXTm = cellfun(@(x) str2num(x),(tempNEXTm(1:end-1,:)));
        end
        NEXTm{i} = tempNEXTm;
        
        
        if isempty(NEXTm{i})
            mSeq{i} = INITm{i};
        else
            mSeq{i} = cat(2,INITm{i},NEXTm{i});
        end
    end
    %%%%%%%%%%%
    
elseif Speller == 2
    
    for i = 1:sum(cell2mat(NumTrials))
        SC{i} = 1:c.NumMatrixColumns.NumericValue;
        SR{i} = (c.NumMatrixColumns.NumericValue+1):2*c.NumMatrixColumns.NumericValue;
    end
    
elseif Speller == 3
    
    for i = 1:sum(cell2mat(NumTrials))
        SC{i} = repmat(7:12,1,NR);
        SR{i} = repmat(1:6,NC,1);
        SR{i} = SR{i}(:)';
    end
    
elseif Speller == 4
    %Changed from this 6/8/17
    %     for i = 1:sum(cell2mat(NumTrials))
    %         SC{i} = repmat(7:12,1,NR);
    %         SR{i} = repmat(1:6,NC,1);
    %         SR{i} = SR{i}(:)';
    %     end
    
    for i = 1:sum(cell2mat(NumTrials))
        SC{i} = 1:size(c.Stimuli.NumericValue,2);
        SR{i} = 1:size(c.Stimuli.NumericValue,2);
    end
    GoalLetterF = c.Stimuli.Value(1,TTSpell(:));
    GoalTrialF = 1:sum(cell2mat(NumTrials));
    ChosenLetterF = SelectedStimulus(SelectedStimulus>0&circshift(SelectedStimulus,[1 0])==0);
end


if Speller == 1
    %Check that the markers saved in the log file are the same as those in the
    %data file
    tStimC = double(nonzeros(StimulusCodeUNArt)); tStimC = tStimC(1:StateDuration:end);
    clear StimCode
    k=1;
    scind = 0;
    for r = 1:length(Runind)
        for i = 1:NumTrials{r}
            clear sctemp
            for j = 1:NumSequences{r}
                %             sctemp{j} = num2cell(tStimC(NumStimCodes*NumSequences{r}*(i-1)+NumStimCodes*(j-1)+1:...
                %                 NumStimCodes*NumSequences{r}*(i-1)+NumStimCodes*(j)));
                sctemp{j} = num2cell(tStimC([1:NumStimCodes]+scind));
                scind = scind+NumStimCodes;
            end
            StimCode{k} = cell2mat(cat(2,sctemp{:}));
            k=k+1;
        end
    end
    
    %StimCode should contain the same values as mSeq
    k = 1;
    clear isgood
    for r = 1:length(Runind)
        for i = 1:NumTrials{r}
            isgood(k) = isequal(StimCode{k},mSeq{k});
            k=k+1;
        end
    end
    if sum(isgood) ~= sum(cell2mat(NumTrials))
        disp('SOMETHING IS WRONG')
    end
end


%If removing data(probably wont do) need to also remove the same portion of
%StimCode

% NumSequences = NumSequences-SequencestoRemove;
% StimCode = StimCode(1:NumSequences,:);

disp(['Done with myfile ' num2str(toc)]);
waitbar(.1,wtbrr);


%% Do CSP
if DoCSP == 1
    if BuildClassifier == 1
        Data_C = (RefFilt*ArtWeights2*Data)';
        [W] = CSP_P300(Data_C,StimulusCodeUNArt,StimulusType,...
            NumChans, NumTrials, NumStimCodes,range,StateDuration);
        [W2] = CSP_P300_ver2(Data_C,SStimulusCodeUNArt,StimulusType,...
            NumChans, NumTrials, NumStimCodes,range,StateDuration);
        disp(['Done with CSP ' num2str(toc)]);
    elseif BuildClassifier == 0
    end
elseif DoCSP == 0
    if BuildClassifier == 1
        W = eye(size(RefFilt,1));
        W2 = eye(size(RefFilt,1));
    else
    end
    disp('CSP not done');
end



%% Visualize Transformed Data
if Figures_On == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Show the time series and spectra of raw,artifacted,filtered,CSPd data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if BuildClassifier == 1
        DataA = ArtWeights2*Data;
        DataAF = RefFilt*DataA;
        DataAFC = W*DataAF;
        % DataAFCT = W'*DataAF;
        DataAFCT = W2*DataAF;
        % DataAFCT = Wother'*DataAF;
        
        
        % HowManyPlots = ceil(NumChans/4);
        % for pp = 1:HowManyPlots
        %     k = 1;
        %     fchan = (pp-1)*4+1:pp*4;
        %     fchane = fchan(fchan<=NumChans);
        %     for ch = fchane
        %         %Check the raw data
        %         figure(pp+3)
        %
        %         rr(k) = subplot(2,2,k); plot(Data(ch,:));
        %         hold on; plot(DataA(ch,:),'r');
        %         plot(DataAF(ch,:),'m');
        %         plot(DataAFCT(ch,:),'g');
        %         plot(StimulusCode,'k');
        %
        %         [pD fD] = pwelch(Data(ch,:),fs*2,fs,fs*2,fs);
        %         [pDA ~] = pwelch(DataA(ch,:),fs*2,fs,fs*2,fs);
        %         [pDAF ~] = pwelch(DataAF(ch,:),fs*2,fs,fs*2,fs);
        %         [pDAFCT ~] = pwelch(DataAFCT(ch,:),fs*2,fs,fs*2,fs);
        %
        %         %Check the spectral content
        %         figure(pp+3+HowManyPlots)
        %         ss(k) = subplot(2,2,k); plot(fD,10*log10(pD),'b'); hold on
        %         plot(fD,10*log10(pDA),'r');
        %         plot(fD,10*log10(pDAF),'m');
        %         plot(fD,10*log10(pDAFCT),'g');
        %
        %         k = k+1;
        %     end
        %     linkaxes(rr,'x');
        %     linkaxes(ss);
        %     set(pp+3,'Name',['Channels ' num2str(fchane)])
        %     set(pp+3+HowManyPlots,'Name',['Spectra, Channels ' num2str(fchane)])
        % end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Plot the P300s for the raw,artifacted,filtered,CSPd data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if FreeSpell == 0
            ClocsOrig = find(StimulusCodeUNArt>0 & StimulusTypeUNArt == 1);
            ClocsOrig = ClocsOrig(1:StateDuration:end);
            if Speller == 4
                NlocsOrig = find(StimulusCodeUNArt>0 & StimulusCodeUNArt<3 & StimulusTypeUNArt == 0);
            else
                NlocsOrig = find(StimulusCodeUNArt>0 & StimulusTypeUNArt == 0);
            end
            NlocsOrig = NlocsOrig(1:StateDuration:end);
            Clocs = find(StimulusCode>0 & StimulusType == 1);
            Clocs = Clocs(1:StateDuration:end);
            if Speller == 4
                Nlocs = find(StimulusCode>0 & StimulusCode<3 & StimulusType == 0);
            else
                Nlocs = find(StimulusCode>0 & StimulusType == 0);
            end
            Nlocs = Nlocs(1:StateDuration:end);
            CCdata = zeros(length(range),size(Data,1),length(ClocsOrig));
            NNdata = zeros(length(range),size(Data,1),length(NlocsOrig));
            CCdataA = zeros(length(range),size(DataA,1),length(Clocs));
            NNdataA = zeros(length(range),size(DataA,1),length(Nlocs));
            CCdataAF = zeros(length(range),size(DataAF,1),length(Clocs));
            NNdataAF = zeros(length(range),size(DataAF,1),length(Nlocs));
            CCdataAFC = zeros(length(range),size(DataAFC,1),length(Clocs));
            NNdataAFC = zeros(length(range),size(DataAFC,1),length(Nlocs));
            CCdataAFCT = zeros(length(range),size(DataAFCT,1),length(Clocs));
            NNdataAFCT = zeros(length(range),size(DataAFCT,1),length(Nlocs));
            toc
            for cl = 1:length(ClocsOrig)
                CCdata(:,:,cl) = Data(:,ClocsOrig(cl)+range-1)';
                if cl<=length(Clocs)
                    CCdataA(:,:,cl) = DataA(:,Clocs(cl)+range-1)';
                    CCdataAF(:,:,cl) = DataAF(:,Clocs(cl)+range-1)';
                    CCdataAFC(:,:,cl) = DataAFC(:,Clocs(cl)+range-1)';
                    CCdataAFCT(:,:,cl) = DataAFCT(:,Clocs(cl)+range-1)';
                end
            end
            toc
            for nl = 1:length(NlocsOrig)
                NNdata(:,:,nl) = Data(:,NlocsOrig(nl)+range-1)';
                if nl<=length(Nlocs)
                    NNdataA(:,:,nl) = DataA(:,Nlocs(nl)+range-1)';
                    NNdataAF(:,:,nl) = DataAF(:,Nlocs(nl)+range-1)';
                    NNdataAFC(:,:,nl) = DataAFC(:,Nlocs(nl)+range-1)';
                    NNdataAFCT(:,:,nl) = DataAFCT(:,Nlocs(nl)+range-1)';
                end
            end
            toc
            
            %             figure('Name','Original Data') %Plot averages
            %             for ch = 1:NumChans
            %                 subplot(4,5,ch); errorbar(mean(CCdata(:,ch,:),3),std(CCdata(:,ch,:),[],3)/sqrt(size(CCdata,3)))
            %                 hold on;
            %                 errorbar(mean(NNdata(:,ch,:),3),std(NNdata(:,ch,:),[],3)/sqrt(size(NNdata,3)),'r')
            %             end
            %
            %             figure('Name','Artifacted Data') %Plot averages
            %             for ch = 1:NumChans
            %                 subplot(4,5,ch); errorbar(mean(CCdataA(:,ch,:),3),std(CCdataA(:,ch,:),[],3)/sqrt(size(CCdataA,3)))
            %                 hold on;
            %                 errorbar(mean(NNdataA(:,ch,:),3),std(NNdataA(:,ch,:),[],3)/sqrt(size(NNdataA,3)),'r')
            %             end
            %
            %             figure('Name','Artifacted and filtered Data') %Plot averages
            %             for ch = 1:NumChans
            %                 subplot(4,5,ch); errorbar(mean(CCdataAF(:,ch,:),3),std(CCdataAF(:,ch,:),[],3)/sqrt(size(CCdataAF,3)))
            %                 hold on;
            %                 errorbar(mean(NNdataAF(:,ch,:),3),std(NNdataAF(:,ch,:),[],3)/sqrt(size(NNdataAF,3)),'r')
            %             end
            %
            %             figure('Name','Artifacted and filtered and CSPed Data') %Plot averages
            %             for ch = 1:NumChans
            %                 subplot(4,5,ch); errorbar(mean(CCdataAFC(:,ch,:),3),std(CCdataAFC(:,ch,:),[],3)/sqrt(size(CCdataAFC,3)))
            %                 hold on;
            %                 errorbar(mean(NNdataAFC(:,ch,:),3),std(NNdataAFC(:,ch,:),[],3)/sqrt(size(NNdataAFC,3)),'r')
            %             end
            %
            %             figure('Name','Artifacted and filtered and CSP''ed Data') %Plot averages
            %             for ch = 1:NumChans
            %                 subplot(4,5,ch); errorbar(mean(CCdataAFCT(:,ch,:),3),std(CCdataAFCT(:,ch,:),[],3)/sqrt(size(CCdataAFCT,3)))
            %                 hold on;
            %                 errorbar(mean(NNdataAFCT(:,ch,:),3),std(NNdataAFCT(:,ch,:),[],3)/sqrt(size(NNdataAFCT,3)),'r')
            %             end
            
            p3p = {CCdata, NNdata};
            p3pn = {'Original Data'};
            if Art ~= 0
                p3p = [p3p; {CCdataA, NNdataA}];
                p3pn = [p3pn, 'Artifacted Data'];
            end
            if RerefVal ~= 0
                p3p = [p3p; {CCdataAF, NNdataAF}];
                p3pn = [p3pn, 'Artifaced and Filtered Data'];
            end
            if DoCSP == 1
                p3p = [p3p; {CCdataAFC, NNdataAFC}; {CCdataAFCT, NNdataAFCT}];
                p3pn = [p3pn, 'Artifacted and filtered and CSPed Data',...
                    'Artifacted and filtered and CSP''ed Data'];
            end
            clrs = [1 0 0;0 0 1;0 1 0;1 0 1];
            if Figures_On == 1
                figure('Name','ERP Of original, artifacted, filtered, CSP data')
                lgnd = {};
                for np = 1:size(p3p,1)
                    cdta = p3p{np,1}; ndta = p3p{np,2};
                    for ch = 1:min([8 NumChans])
                        subplot(2,4,ch);
                        %                         plot(range/fs,squeeze(cdta(:,ch,:)),'Color',[.7 .7 1]);
                        %                         hold on
                        %                         plot(range/fs,squeeze(ndta(:,ch,:)),'Color',[1 .7 .7]);
                        plot(range/fs,nanmean(cdta(:,ch,:),3),'Color',clrs(np,:));
                        hold on
                        plot(range/fs,nanmean(ndta(:,ch,:),3),'Color',.5*clrs(np,:));
                    end
                    lgnd = [lgnd [p3pn{np} '_{Choice}'] [p3pn{np} '_{Not}']];
                end
                legend(lgnd)
            end
            toc
            
            
            % %For NIH Grant-used no reref, no art, yes csp -- plotted by running through
            % %each subject, plotting, saving, and adding to the previous plot.
            % %saved as NIHgrantP300_st_an_sm_aj.fig
            % %S1 - AndrewP3_14001 R03/04 -- csp chan 1
            % %S2 - SumithraP3001 R03/04 -- csp chan 1
            % %S3 - AnjumGridCB001 R01 -- csp chan 1
            % %S4 - SteveGridCB002 R03 -- csp chan 1
            % ch=1;
            % figure
            % hold on
            % plot(range/fs,mean(CCdataAFC(:,ch,:),3)); hold on;
            % plot(range/fs,mean(NNdataAFC(:,ch,:),3),'r'); xlim([0 range(end)/fs])
            % xlabel('Time (s)'); ylabel('Amplitude \muV');
            
            
            
            
        end
    end
    
    clear DataA DataAF DataAFC DataAFCT CCData CCdataA CCdataAF CCdataAFC
    clear CCdataAFCT NNdata NNdataA NNdataAF NNdataAFC NNdataAFCT
    
    
    %StimulusType are repeated
    %
    % [aa bb] = butter(5,[.5 30]/128);
    % aaa = filtfilt(aa,bb,Data);
    %
    % figure
    % plot(Data(:,1)); hold on; plot(aaa(:,1),'r');
    %
    % figure
    % pwelch(Data(:,1),1024,0,512,256); hold on;
    % pwelch(aaa(:,1),1024,0,512,256);
    %
    % figure
    % plot(StimulusCode);
    % hold on;
    % plot(StimulusType,'r');
    
    
    %     %Reshape data so it can be multiplied by Wp
    %     Data2=reshape(AvgEEGData,NumStimCodes*NumTrials,NumChans,size(AvgEEGData,2));
    %     Z = zeros(size(Data2,1),size(Wp,1),size(Data2,3));
    %     %Transform data with spatial filter Wp
    %     for tt = 1:size(Data2,1)
    %         Z(tt,:,:) = Wp*squeeze(Data2(tt,:,:));
    %     end
    %     %Reshape data
    %     Z2 = reshape(Z,NumStimCodes*NumTrials*NumChans,size(Data2,3));
    
    
    % P300 sanity check (2)
end

disp(['Done with vis ' num2str(toc)]);
waitbar(.15,wtbrr);

%% Save/Load Spatial Filter
if BuildClassifier == 1
    
    %Which CSP channels to Keep?
    
    
    %     ChansToKeep = input('Which Channels to keep?'
    if DoCSP == 1
        if RerefVal == 1 %This is because the
            %last CSP channel unser CAR rereferencing is zeros.
            if NumChans <= 4
                ChansToKeep = 1:NumChans-1;
            else
                ChansToKeep = [1 2 NumChans-2 NumChans-1];
            end
        else
            if NumChans <= 3
                ChansToKeep = 1:NumChans;
            else
                ChansToKeep = [1 2 NumChans-1 NumChans];
            end
        end
        
    else
        ChansToKeep = 1:size(W,1);
    end
    
    NumChans = length(ChansToKeep);
    W = W(ChansToKeep,:);
    
    
    Data(isnan(Data))=0;
    SpatFilt = W*RefFilt*ArtWeights2;
    %     SpatFilt2 = W2*RefFilt*ArtWeights2;
    %Combine Spatial Filter and Artifact Matrix into Spatial Filter for online
    %use.  Here i am first doing artifact rejection, then spatial
    %filtering.  There is a slight difference with the order.
    
    %Save the Spatial Filter
    %     dlmwrite(['C:\Documents and Settings\amg5106\My Documents\'...
    %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' Name Session '\'...
    %         Name 'S' Session 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
    try
        dlmwrite([bcifold '\data\' Name Sessions{max_sess_i} '\'...
            Name 'S' Sessions{max_sess_i} 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
    catch
        dlmwrite(['P:\ALS Proj Data\' Name Sessions{max_sess_i} '\'...
            Name 'S' Sessions{max_sess_i} 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
    end
    DD = SpatFilt*Data;
elseif BuildClassifier == 0
    try
        SpatFilt = dlmread([NameBase 'SpatialFilter.txt'],'\t');
    catch
        SpatFilt = dlmread([NameBase2 'SpatialFilter.txt'],'\t');
    end
    
    
    %     disp('DONT FORGET TO CHANGE THIS BACK')
    %     SpatFilt = eye(size(Data,1));
    
    Data(isnan(Data))=0;
    DD = SpatFilt*Data;
    NumChans = size(DD,1);
    
    
    
    
end

disp(['Done with filtering ' num2str(toc)]);
waitbar(.2,wtbrr);
%% EyeTracker
if Speller==4 %Can include this later -- for now just ignore eye tracking for audio speller
    eyye = [];
else
    if sum(~isnan(GazeX))~=0
        [gazeacc, gazevar, gazeinvalid] = P300_EyeTracking_inP300Classifier_v3(...
            GazeX,GazeY, Data(EOGloc,:), SC,SR, StimulusCodeUNArt, ...
            SequencePhase, StimulusTypeUNArt, StateDuration, ...
            NumTrials,NumSequences, Speller,...
            fs, c, Figures_On, 0, TTSpell, AppPos);
        disp(['Finished with eye tracking ' num2str(toc)]);
    else
        disp('No eye tracking data detected');
    end
    
    waitbar(.25,wtbrr);
    eyye.gazeacc = gazeacc;
    eyye.gazevar = gazevar;
    eyye.gazeinvalid = gazeinvalid;
end
%% Process and Organize Data



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Organize Data into choice and not-choice targets


if FreeSpell==0
    targetstim = double(StimulusCodeUNArt(logical(StimulusTypeUNArt)));
    targetstim = targetstim(1:StateDuration:length(targetstim));
    
    k = 1;
    ki = 1;
    for r = 1:length(NumTrials)
        for i = 1:NumTrials{r}
            if Speller == 4
                ts(ki,:) = unique(targetstim(k:k+NumSequences{r}-1))';
                k=k+NumSequences{r};
            else
                ts(ki,:) = unique(targetstim(k:k+NumSequences{r}*2-1))';
                k=k+NumSequences{r}*2;
            end
            ki=ki+1;
        end
    end
    targetstim = ts;
    
else
    targetstim = repmat([1 2],sum(cell2mat(NumTrials)),1); %Put fake target stims in... just
    %so the code below runs and we group the data for classification.
    %We will ignore this when classifying later
end


%%%% 2/3/17 IS THIS EVEN USEFUL ANYMORE????
%Initialize counters for each stimulus code
Tind = ones(1,NumStimCodes);
maxseq = max(cell2mat(NumSequences));
AllData = zeros(NumStimCodes*(sum(cell2mat(NumTrials))*maxseq)*NumChans,length(range));
AllData_c = zeros(NumStimCodes*(sum(cell2mat(NumTrials))*maxseq)*NumChans,1);
AllData_t = zeros(NumStimCodes*(sum(cell2mat(NumTrials))*maxseq)*NumChans,1);
ad_i = 1;
t_ind = 1;
maxseq = max(cell2mat(NumSequences));
seqdiff = maxseq-cell2mat(NumSequences);
TrialTargData = cell(sum(cell2mat(NumTrials)),NumStimCodes);
for r = 1:length(NumTrials)
    
    tic
    for i = 1:NumTrials{r}
        
        for j = 1:NumStimCodes %Loop through stimulus Codes
            
            TargLoc = find(StimulusCode == j);
            if t_ind == sum(cell2mat(NumTrials))
                TargLoc(TargLoc<TrialStartTime(t_ind)) = [];
            else
                TargLoc(TargLoc<TrialStartTime(t_ind)|TargLoc>TrialStartTime(t_ind+1)) = [];
            end
            TrialTargLoc{t_ind,j} = TargLoc(1:StateDuration:end);
            
            if ismember(j,targetstim(i,:))
                tmrk = 1;
            else
                tmrk = 0;
            end
            %         %Now take the locations of the first NumSeq of each of these vectors, and
            %         %increment a counter by NumSeq
            %         TrialTargLoc{i,j} = TargLoc(Tind(j):Tind(j)+NumSequences-1);
            %         Tind(j) = Tind(j)+NumSequences+SequencestoRemove;
            %
            for ch = 1:NumChans %Loop through channels  %CAN PROBABLY REMOVE THIS LOOP
                TrialTargData{t_ind,j}{ch} = [];
                TrialTargData_i{t_ind,j}{ch} = [];
                %             figure
                
                
                sst = length(TrialTargLoc{t_ind,j});
                TrialTargData{t_ind,j}{ch} = reshape([DD(ch,repmat(TrialTargLoc{t_ind,j},1,length(range))'+repmat(range-1,sst,1)') ...
                    NaN(1,(maxseq-sst)*length(range))],length(range),maxseq)';
                TrialTargData_i{t_ind,j}{ch} = (ad_i:ad_i+maxseq-1)';
                AllData(ad_i:ad_i+maxseq-1,:) = TrialTargData{t_ind,j}{ch};
                AllData_c(ad_i:ad_i+maxseq-1,:) = ch;
                AllData_t(ad_i:ad_i+maxseq-1,:) = tmrk;
                ad_i = ad_i+maxseq;
                
                %                 s_i = 0;
                %                 for s = 1:length(TrialTargLoc{t_ind,j})%MAY BE ABLE TO REMOVE THIS LOOP %Loop through sequences per trial, extract data 0-.6 seconds after stimulus
                %
                %                     TrialTargData{t_ind,j}{ch} = [TrialTargData{t_ind,j}{ch}; DD(ch,TrialTargLoc{t_ind,j}(s)+range-1)];
                %                     TrialTargData_i{t_ind,j}{ch} = [TrialTargData_i{t_ind,j}{ch}; ad_i];
                %                     AllData(ad_i,:) = DD(ch,TrialTargLoc{t_ind,j}(s)+range-1);
                %                     AllData_c(ad_i) = ch;
                %                     AllData_t(ad_i) = tmrk;
                %                     ad_i = ad_i+1;
                %                     s_i = s_i+1;
                %                 end
                %
                % %                 Add in NaNs for remaining sequences less than the maximum number of
                % %                 sequences, just so all sessions have the same number
                %                 TrialTargData{t_ind,j}{ch} = [TrialTargData{t_ind,j}{ch}; NaN(maxseq-s_i,length(range))];
                %                 TrialTargData_i{t_ind,j}{ch} = [TrialTargData_i{t_ind,j}{ch}; (ad_i:ad_i+(maxseq-s_i)-1)'];
                %                 AllData(ad_i:ad_i+(maxseq-s_i)-1,:) = NaN(maxseq-s_i,length(range));
                %                  AllData_c(ad_i:ad_i+(maxseq-s_i)-1,:) = ch;
                %                  AllData_t(ad_i:ad_i+(maxseq-s_i)-1,:) = tmrk;
                %                  ad_i = ad_i+(maxseq-s_i);
            end
        end
        t_ind=t_ind+1;
    end
    toc
end
clear DD Data
disp(['Time for organizing data ' num2str(toc)]);
waitbar(.3,wtbrr);


%Find indices of bad (unaveraged) epochs here, then average sequences, for
%each channel separately and for choice and non-choice separately.
BadEpoch = zeros(NumStimCodes*(double(cell2mat(NumTrials))*double(cell2mat(NumSequences))')*NumChans,1);
if BuildClassifier == 1
    for ch = 1:NumChans
        
        tempind_c = find(AllData_c == ch & AllData_t == 1);
        tempdata_c = AllData(tempind_c,:);
        tempm_c = mean(tempdata_c); temps_c = std(tempdata_c);
        overtemp_c = sum(tempdata_c > repmat(tempm_c+5*temps_c,size(tempdata_c,1),1),2);
        undertemp_c = sum(tempdata_c < repmat(tempm_c-5*temps_c,size(tempdata_c,1),1),2);
        badind_c = (overtemp_c+undertemp_c)>0;
        BadEpoch(tempind_c) = badind_c;
        
        tempind_n = find(AllData_c == ch & AllData_t == 0);
        tempdata_n = AllData(tempind_n,:);
        tempm_n = mean(tempdata_n); temps_n = std(tempdata_n);
        overtemp_n = sum(tempdata_n > repmat(tempm_n+5*temps_n,size(tempdata_n,1),1),2);
        undertemp_n = sum(tempdata_n < repmat(tempm_n-5*temps_n,size(tempdata_n,1),1),2);
        badind_n = (overtemp_n+undertemp_n)>0;
        BadEpoch(tempind_n) = badind_n;
    end
end
clear AllData

%%v3
% %Find indices of bad (unaveraged) epochs here, then average sequences, for
% %each channel separately and for choice and non-choice separately.
% BadEpoch = zeros(NumStimCodes*(double(cell2mat(NumTrials))*double(cell2mat(NumSequences))')*NumChans,1);
% if BuildClassifier == 1
%     for ch = 1:NumChans
%
%         tempind_c = find(AllData_c == ch & AllData_t == 1);
%         tempdata_c = AllData(tempind_c,:);
%         tempm_c = mean(tempdata_c); temps_c = std(tempdata_c);
%         overtemp_c = sum(tempdata_c > repmat(tempm_c+5*temps_c,size(tempdata_c,1),1),2);
%         undertemp_c = sum(tempdata_c < repmat(tempm_c-5*temps_c,size(tempdata_c,1),1),2);
%         badind_c = (overtemp_c+undertemp_c)>0;
%         BadEpoch(tempind_c) = badind_c;
%
%         tempind_n = find(AllData_c == ch & AllData_t == 0);
%         tempdata_n = AllData(tempind_n,:);
%         tempm_n = mean(tempdata_n); temps_n = std(tempdata_n);
%         overtemp_n = sum(tempdata_n > repmat(tempm_n+5*temps_n,size(tempdata_n,1),1),2);
%         undertemp_n = sum(tempdata_n < repmat(tempm_n-5*temps_n,size(tempdata_n,1),1),2);
%         badind_n = (overtemp_n+undertemp_n)>0;
%         BadEpoch(tempind_n) = badind_n;
%     end
% end
% clear AllData

NumSequences = maxseq;
NumTrials = sum(cell2mat(NumTrials));
NumStimCodes = double(NumStimCodes);
EEGData = zeros(NumStimCodes*NumTrials*NumSequences*NumChans,length(range));
EEGDataRR = zeros(NumStimCodes*NumTrials*NumSequences*NumChans,length(range));
% EEGData = zeros(NumStimCodes*(double(cell2mat(NumTrials))*double(cell2mat(NumSequences))')*NumChans,length(range));
% EEGDataRR = zeros(NumStimCodes*(double(cell2mat(NumTrials))*double(cell2mat(NumSequences))')*NumChans,length(range));
ChanLabel = zeros(size(EEGData,1),1);
ClassLabel = zeros(size(EEGData,1),1);
AvgChanLabel = zeros(NumTrials*NumStimCodes*NumChans,1);
AvgClassLabel = zeros(NumTrials*NumStimCodes*NumChans,1);
idx = cell(NumChans,NumTrials,NumStimCodes);
idx = cellfun(@(x) 0, idx, 'UniformOutput', false);
for  ns = 1:NumSequences
    AvgEEGData{ns} = zeros(NumTrials*NumStimCodes*NumChans,length(range));
    AvgEEGDataRR{ns} = zeros(NumTrials*NumStimCodes*NumChans,length(range));
    for ch = 1:NumChans %Loop through channels
        for i = 1:NumTrials %Loop through Trials
            for j = 1:NumStimCodes
                
                
                trl_index = TrialTargData_i{i,j}{ch};
                good_trl = find(BadEpoch(trl_index)==0);
                
                
                if ns == NumSequences
                    idx = (ch-1)*(NumTrials*NumStimCodes*NumSequences)+(i-1)*...
                        (NumStimCodes*NumSequences)+(j-1)*NumSequences+1:(ch-1)*...
                        (NumTrials*NumStimCodes*NumSequences)+(i-1)*(NumStimCodes*...
                        NumSequences)+j*NumSequences;
                    EEGData(idx,:) = [TrialTargData{i,j}{ch}(good_trl,:); NaN(NumSequences-length(good_trl),length(range))];
                    EEGDataRR(idx,:) =[TrialTargData{i,j}{ch}; NaN(NumSequences-length(trl_index),length(range))];
                    
                end
                
                if ns>length(good_trl)
                    epc_to_avg = good_trl;
                else
                    epc_to_avg = good_trl(1:ns);
                end
                
                avg_index = (ch-1)*(NumTrials*NumStimCodes)+(i-1)*...
                    (NumStimCodes)+(j-1)+1:(ch-1)*(NumTrials*NumStimCodes)+...
                    (i-1)*(NumStimCodes)+(j);
                
                if isempty(epc_to_avg)
                    AvgEEGDataRR{ns}(avg_index,:) = NaN(1,length(range));
                    AvgEEGDataRR{ns}(avg_index,:) = NaN(1,length(range));
                else
                    AvgEEGData{ns}(avg_index,:) = nanmean(TrialTargData{i,j}{ch}(epc_to_avg,:),1);
                    AvgEEGDataRR{ns}(avg_index,:) = nanmean(TrialTargData{i,j}{ch});
                end
                %             if ch==1 && i==1 && j==1
                %             disp('Detrending P300 epochs!');
                %             end
                %                              EEGData(index,:) = detrend(TrialTargData{i,j}{ch},0);
                %                              AvgEEGData(avg_index,:) = detrend(mean(TrialTargData{i,j}{ch},1),0);
                %
                
                %             SScodes(index)=j;
                if ns == NumSequences
                    if ismember(j,targetstim(i,:)) %Choice targets = 1
                        ChanLabel(idx) = repmat(ch,NumSequences,1);
                        ClassLabel(idx) = ones(NumSequences,1);
                        AvgChanLabel(avg_index) = ch;
                        AvgClassLabel(avg_index) = 1;
                    else
                        if Speller == 4 && j>2 %Class beeps in the audio speller separately
                            ChanLabel(idx) = repmat(ch,NumSequences,1);
                            ClassLabel(idx) = 3*ones(NumSequences,1);
                            AvgChanLabel(avg_index) = ch;
                            AvgClassLabel(avg_index) = 3;
                        else
                            ChanLabel(idx) = repmat(ch,NumSequences,1);
                            ClassLabel(idx) = 2*ones(NumSequences,1);
                            AvgChanLabel(avg_index) = ch;
                            AvgClassLabel(avg_index) = 2;
                        end
                    end
                end
                
                
                %             dataloc{ch}(index) = TrialTargLoc{i,j};
            end
            
        end
    end
end

if FreeSpell == 1
    clear ClassLabel AvgClassLabel
end
disp(['Time to rearrange data (2) ' num2str(toc)]);
waitbar(.35,wtbrr);
if Figures_On == 1    figure
    tt(1) = subplot(121);imagesc(EEGData); colorbar
    tt(2) = subplot(122); imagesc(EEGDataRR); colorbar;
    linkaxes(tt)
    
end


%Plot averages before and after epoch removal
tic
if Figures_On == 1
    if Avv == 1
        if FreeSpell == 0
            figure('Name','Target and non-target ERPs of filtered data, before and after epoch removal')
            for ch = 1:NumChans
                subplot(4,5,ch);
                plot(range/fs,mean(AvgEEGDataRR{ns}(AvgClassLabel==2&AvgChanLabel==ch,:)),'Color',[1 .7 .7]);
                hold on;
                plot(range/fs,mean(AvgEEGDataRR{ns}(AvgClassLabel==1&AvgChanLabel==ch,:)),'Color',[.7 .7 1]);
                if Speller == 4
                    plot(range/fs,mean(AvgEEGDataRR{ns}(AvgClassLabel==3&AvgChanLabel==ch,:)),'Color',[.7 .7 .7]);
                end
                plot(range/fs,mean(AvgEEGData{ns}(AvgClassLabel==2&AvgChanLabel==ch,:)),'Color','r');
                plot(range/fs,mean(AvgEEGData{ns}(AvgClassLabel==1&AvgChanLabel==ch,:)),'Color','b');
                if Speller == 4
                    plot(range/fs,mean(AvgEEGData{ns}(AvgClassLabel==3&AvgChanLabel==ch,:)),'Color','k');
                end
                title(ChannelNames{ch});
            end
        elseif FreeSpell == 1
            figure('Name','Averaged ERPs of filtered data')
            for ch = 1:NumChans
                subplot(4,5,ch);
                plot(range,AvgEEGDataRR{ns}(AvgChanLabel==ch,:));
            end
        end
    else %Avv=0
        if FreeSpell == 0
            figure('Name','Target and non-target ERPs of filtered data, before and after epoch removal')
            
            t1=sum(isnan(EEGData),2)>0;
            t2=sum(isnan(EEGDataRR),2)>0;
            badremoved = zeros(size(EEGDataRR));
            badremoved(t1,:) = EEGDataRR(t1,:);
            
            
            for ch = 1:NumChans
                subplot(4,5,ch);
                plot(range/fs,nanmean(EEGDataRR(ClassLabel==2&ChanLabel==ch,:)),'Color',[1 .7 .7]);
                hold on;
                plot(range/fs,nanmean(EEGDataRR(ClassLabel==1&ChanLabel==ch,:)),'Color',[.7 .7 1]);
                plot(range/fs,nanmean(badremoved(ClassLabel==2&ChanLabel==ch,:)),'Color',[.7 1 .7]);
                if Speller == 4
                    plot(range/fs,nanmean(badremoved(ClassLabel==3&ChanLabel==ch,:)),'Color',[.7 .7 .7]);
                end
                plot(range/fs,nanmean(EEGData(ClassLabel==2&ChanLabel==ch,:)),'Color','r');
                plot(range/fs,nanmean(EEGData(ClassLabel==1&ChanLabel==ch,:)),'Color','b');
                if Speller == 4
                    plot(range/fs,nanmean(EEGData(ClassLabel==3&ChanLabel==ch,:)),'Color','k');
                end
                plot(range/fs,nanmean(badremoved(ClassLabel==1&ChanLabel==ch,:)),'Color','g');
                
                title(ChannelNames{ch});
            end
        elseif FreeSpell == 1
            figure('Name','Averaged ERPs of filtered data')
            for ch = 1:NumChans
                subplot(4,5,ch);
                plot(range,nanmean(EEGDataRR(ChanLabel==ch,:)),'b');
            end
        end
    end
end

clear EEGDataRR AvgEEGDataRR


disp(['Plot after epoch removal ' num2str(toc)]);

%% Data Reduction - for LDA Classifiers
if ReducT == 1
    %Decimate Data - use decimated data in classifier
    %
    %     clear dEEGData dAvgEEGData
    %     for dd = 1:size(EEGData,1)
    %         dEEGData(dd,:) = decimate(EEGData(dd,range),8);
    %     end
    %     for dd = 1:size(AvgEEGData,1)
    %         dAvgEEGData(dd,:) = decimate(AvgEEGData(dd,range),8);
    %     end
    %     LDAtimes = round(decimate(range,8));
    %
    %
    %     inpch = [];
    %     if Avv == 0
    %         TrainD = [];
    %         for ch = 1:max(ChanLabel)
    %             TrainD = [TrainD double(dEEGData(ChanLabel==ch,:))];
    %             inpch = [inpch ch*ones(1,size(dEEGData,2))];
    %         end
    %     else
    %         AvgTrainD = [];
    %         for ch=1:max(AvgChanLabel)
    %             AvgTrainD = [AvgTrainD double(dAvgEEGData(AvgChanLabel==ch,:))];
    %             inpch = [inpch ch*ones(1,size(dAvgEEGData,2))];
    %         end
    %     end
    
elseif ReducT == 2
    %     %Find timepoints where averages are significantly different - pick from
    %     %these to be the features of the LDA.
    %     for ch = 1:max(ChanLabel)
    %         %             figure
    %         Cch = ChanLabel==ch;
    %         Cn = ClassLabel=='N';
    %         Cc = ClassLabel=='C';
    %         nloc = Cch&Cn;
    %         cloc = Cch&Cc;
    %         %             plot(mean(EEGData(nloc,range)),'r');
    %         %             hold on;
    %         %             plot(mean(EEGData(cloc,range)),'b');
    %         for i = 1:length(range)
    %             [tstat(ch,i) pval(ch,i)] =ttest(EEGData(cloc,range(i)),mean(EEGData(nloc,range(i))));
    %         end
    %         [~, LDAtimes(ch)]=min(pval(ch,:));
    %         dEEGData(Cch,:) = EEGData(Cch,LDAtimes(ch));
    %
    %
    %         %             figure
    %         Cch = AvgChanLabel==ch;
    %         Cn = AvgClassLabel=='N';
    %         Cc = AvgClassLabel=='C';
    %         nloc = Cch&Cn;
    %         cloc = Cch&Cc;
    %         %             plot(mean(AvgEEGData(nloc,range)),'r');
    %         %             hold on;
    %         %             plot(mean(AvgEEGData(cloc,range)),'b');
    %         for i = 1:length(range)
    %             [Atstat(ch,i) Apval(ch,i)] =ttest(AvgEEGData(cloc,range(i)),mean(AvgEEGData(nloc,range(i))));
    %         end
    %         [~, AvgLDAtimes(ch)]=min(Apval(ch,:));
    %         dAvgEEGData(Cch,:) = AvgEEGData(Cch,AvgLDAtimes(ch));
    %
    %     end
    %
    %     inpch = [];
    %     if Avv == 0
    %         TrainD = [];
    %         for ch = 1:max(ChanLabel)
    %             TrainD = [TrainD double(dEEGData(ChanLabel==ch,:))];
    %             inpch = [inpch ch*ones(1,size(dEEGData,2))];
    %         end
    %     elseif Avv == 1
    %         AvgTrainD = [];
    %         for ch=1:max(AvgChanLabel)
    %             AvgTrainD = [AvgTrainD double(dAvgEEGData(AvgChanLabel==ch,:))];
    %             inpch = [inpch ch*ones(1,size(dAvgEEGData,2))];
    %         end
    %     end
    
elseif ReducT == 3 %%Low pass filter, downsample to 20Hz
    
    %lowpass filter at 20 Hz
    LPv = 20;
    [bb2 aa2] = butter(5,LPv/(fs/2),'low');
    LDAtimes = downsample(range,fix(fs/(LPv)));
    
    
    inpch = [];
    x_filtered = filtfilt(bb2,aa2,EEGData');
    x_down = downsample(x_filtered,fix(fs/(LPv)))';
    
    TrainD = []; NotDownD = [];
    for ch = 1:max(ChanLabel)
        TrainD = [TrainD double(x_down(ChanLabel==ch,:))];
        inpch = [inpch ch*ones(1,size(x_down,2))];
        %         NotDownD = [NotDownD double(x_filtered(:,ChanLabel==ch)')];
    end
    for ns = 1:NumSequences
        
        
        
        x_filtered = filtfilt(bb2,aa2,AvgEEGData{ns}');
        x_down = downsample(x_filtered,fix(fs/(LPv)))';
        %         x_down2 = decimate(x_filtered,fix(fs/LPv))';
        
        %                 figure
        %                 plot(AvgEEGData(142,:));
        %                 hold on;
        %                 plot(x_filtered(:,142),'g')
        %                 plot(LDAtimes,x_down(142,:),'r');
        AvgTrainD{ns} = [];
        for ch=1:max(AvgChanLabel)
            AvgTrainD{ns} = [AvgTrainD{ns} double(x_down(AvgChanLabel==ch,:))];
        end
    end
    
    
    
elseif ReducT == 0
    LDAtimes = range;
    for ns = 1:NumSequences
        inpch = [];
        
        TrainD = [];
        for ch = 1:max(ChanLabel)
            TrainD = [TrainD double(EEGData(ChanLabel==ch,:))];
            inpch = [inpch ch*ones(1,size(EEGData,2))];
        end
        
        AvgTrainD{ns} = [];
        for ch=1:max(AvgChanLabel)
            AvgTrainD{ns} = [AvgTrainD{ns} double(AvgEEGData{ns}(AvgChanLabel==ch,:))];
        end
        
    end
    
    
end

disp(['Time to data reduce ' num2str(toc)]);
waitbar(.4,wtbrr);
%% P300 Sanity Check (3)
if Figures_On == 1
    if FreeSpell == 0
        figure('Name','Average target and non-target ERPs of filtered, down sampled, and epoch rejected data');
        for ch = 1:NumChans
            subplot(4,5,ch);
            plot(LDAtimes,AvgTrainD{ns}(AvgClassLabel(AvgChanLabel==1)==2,...
                (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes))','Color',[1 .7 .7]); hold on;
            plot(LDAtimes,AvgTrainD{ns}(AvgClassLabel(AvgChanLabel==1)==1,...
                (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes))','Color',[.7 .7 1]);
            if Speller == 4
                plot(LDAtimes,AvgTrainD{ns}(AvgClassLabel(AvgChanLabel==1)==3,...
                    (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes))','Color',[.7 .7 .7]);
            end
            plot(LDAtimes,nanmean(AvgTrainD{ns}(AvgClassLabel(AvgChanLabel==1)==2,...
                (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes))),'r'); hold on;
            plot(LDAtimes,nanmean(AvgTrainD{ns}(AvgClassLabel(AvgChanLabel==1)==1,...
                (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes))),'b');
            if Speller == 4
                plot(LDAtimes,nanmean(AvgTrainD{ns}(AvgClassLabel(AvgChanLabel==1)==3,...
                    (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes))),'k');
            end
        end
    elseif FreeSpell == 1
    end
end
disp(['Time to plot (2) ' num2str(toc)]);

%% Classification
ClassTypes = {'SW','LDA','SS','LO','rLDA','Manual'};
if BuildClassifier == 1
    
    %Four possible parameter files 8/1/16
    %BCI Trainer
    if ctype == 1, ptype = 'bcitraining'; end;
    %Symbolic Speller
    if ctype == 2, ptype = 'symbolic'; end;
    %Dual Speller
    if ctype == 3, ptype = 'dual'; end;
    %Notepad Speller
    if ctype == 4, ptype = 'notepad'; end;
    %Audio Speller
    if ctype == 5, ptype = 'audio'; end;
    
    if ctype == 5
        pfolder = ls([bcifold '\parms\AudioSpell*']);
        pfolder2 = [bcifold '\parms\' pfolder '\audiospeller_' capt '_2monitor.prm'];
    else
        pfolder = ls([bcifold '\parms\CBSpell*']);
        pfolder2 = [bcifold '\parms\' pfolder '\' ptype '_' capt '_40targface_2monitor.prm'];
    end
    fpm = fopen(pfolder2,'r');
    dmmy = 0;
    if fpm == -1
        disp('Parameter file does not exist!');
        fpm = fopen([bcifold  '\parms\' pfolder '\dummyfile.prm'],'r');
        dmmy = 1;
    end
    
    fpm2=textscan(fpm,'%s','delimiter','\n');
    fclose(fpm);
    for ll = 1:length(fpm2{1})
        LinInd(ll) = ~isempty(strfind(fpm2{1}{ll},'Filtering:Linear'));
        SNInd(ll) = ~isempty(strfind(fpm2{1}{ll},'SubjectName'));
        SEInd(ll) = ~isempty(strfind(fpm2{1}{ll},'SubjectSession'));
        SFInd(ll) = ~isempty(strfind(fpm2{1}{ll},'SpatialFilter matrix'));
        AccEvInd(ll) = ~isempty(strfind(fpm2{1}{ll},'AccumulateEvidence'));
        MinEvInd(ll) = ~isempty(strfind(fpm2{1}{ll},'MinimumEvidence'));
        DispRes(ll) = ~isempty(strfind(fpm2{1}{ll},'DisplayResults'));
        if ctype == 5
            TTS(ll) = ~isempty(strfind(fpm2{1}{ll},'ToBeCopied'));
        else
            TTS(ll) = ~isempty(strfind(fpm2{1}{ll},'TextToSpell'));
        end
        ImpInd(ll) = ~isempty(strfind(fpm2{1}{ll},'AcquisitionMode'));
        IntInd(ll) =  ~isempty(strfind(fpm2{1}{ll},'InterpretMode'));
        STVYInd(ll) = ~isempty(strfind(fpm2{1}{ll},'Sensitivity'));
    end
    LinInd = find(LinInd==1);
    SNInd = find(SNInd==1); SNInd = SNInd(1);  %Needed to add for 4903BCI
    SEInd = find(SEInd==1); SEInd = SEInd(1);  %Needed to add for 4903BCI
    SFInd = find(SFInd==1);
    AccEvInd = find(AccEvInd==1);
    MinEvInd = find(MinEvInd==1);
    DispResInd = find(DispRes==1);
    TTSInd = find(TTS==1);
    ImpInd = find(ImpInd==1);
    SBegInd = strfind(fpm2{1}{SNInd},'='); SEndInd = strfind(fpm2{1}{SNInd},'Name %');
    SEBegInd = strfind(fpm2{1}{SEInd},'='); SEEndInd = strfind(fpm2{1}{SEInd},'% %');
    SFBegInd = strfind(fpm2{1}{SFInd},'='); SFEndInd = strfind(fpm2{1}{SFInd},'//');
    BegInd = strfind(fpm2{1}{LinInd},'}'); EndInd = strfind(fpm2{1}{LinInd},'//');
    BegInd2 = strfind(fpm2{1}{LinInd},'= '); EndInd2 = strfind(fpm2{1}{LinInd},' {');
    DispResBegInd = strfind(fpm2{1}{DispResInd},'= '); DispResEndInd = strfind(fpm2{1}{DispResInd},'//');
    TTSBegInd = strfind(fpm2{1}{TTSInd},'= '); TTSEndInd = strfind(fpm2{1}{TTSInd},'//');
    if ~isempty(AccEvInd)
        ABegInd = strfind(fpm2{1}{AccEvInd},'=');
        MBegInd = strfind(fpm2{1}{MinEvInd},'='); MEndInd = strfind(fpm2{1}{MinEvInd},'//');
    end
    if ~isempty(ImpInd)
        ImpBegInd = strfind(fpm2{1}{ImpInd},'= '); ImpEndInd = strfind(fpm2{1}{ImpInd},'//');
    end
    STVYInd = find(STVYInd==1);
    if ~isempty(STVYInd)
        STVYBegInd = strfind(fpm2{1}{STVYInd},'= '); STVYEndInd = strfind(fpm2{1}{STVYInd},'//');
    end
    
    %      NumCrossVal = 10;
    %
    %         for ccc = 1:NumCrossVal
    %             test_ind = (ccc-1)*size(AvgTrainD,1)/NumCrossVal+1:...
    %                 ccc*size(AvgTrainD,1)/NumCrossVal;
    %             train_ind = find(~ismember(1:size(AvgTrainD),test_ind));
    %             cv_train = AvgTrainD(train_ind,:);
    %             cv_test = AvgTrainD(test_ind,:);
    %             cv_cl = AvgClassLabel(train_ind);
    %
    %
    %             if Cfier == 1
    %                 %Stepwise Linear Regression
    %                 if Avv == 0
    %                     [B, SE, PVAL, INMODEL] = stepwisefit(TrainD,ClassLabel(ChanLabel==1));
    %                 else
    %                     [B, SE, PVAL, INMODEL] = stepwisefit(cv_train,cv_cl);
    %                 end
    %                 SigLocs = find(INMODEL==1);
    %                 Classweight = -B;
    %
    %             elseif Cfier == 2
    %                 %LDA
    %                 %Classifier Built from individual trials (dont think this is what I
    %                 %want to do)
    %                 if Avv == 0
    %                     [class,err,post,logp,coeff] = classify(TrainD,TrainD,ClassLabel(ChanLabel==1));
    %                 elseif Avv == 1
    %                     [class,err,post,logp,coeff] = classify(AvgTrainD,AvgTrainD,AvgClassLabel(AvgChanLabel==1));
    %                 end
    %                 Classweight = coeff(1,2).linear;
    %                 SigLocs = find(Classweight~=0);
    %
    %
    %             elseif Cfier == 3
    %                 % Steve's LDA
    %                 %This produces a crazy gamma matrix similar to the regular LDA, which
    %                 %doesnt work.  Figure out how to make regular LDA output a useful gamma
    %                 %matrix and then come back to this
    %
    %                 if Avv == 0
    %                     cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
    %                         'ALS Project\Part2 - P300\P300Classifier\LDA Code\Redo'])
    %                     [Gam,z] = runLDA('Andrew',TrainD,ClassLabel(ChanLabel==1));
    %                     cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
    %                         'ALS Project\Part2 - P300\P300Classifier'])
    %                 elseif Avv == 1
    %                     cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
    %                         'ALS Project\Part2 - P300\P300Classifier\LDA Code\Redo'])
    %                     [Gam,z] = runLDA('Andrew',AvgTrainD,AvgClassLabel(AvgChanLabel==1));
    %                     cd (['C:\Documents and Settings\amg5106\My Documents\SugarSync\Year 4\'...
    %                         'ALS Project\Part2 - P300\P300Classifier'])
    %                 end
    %                 Classweight = Gam;
    %                 SigLocs = find(Gam~=0);
    %             end
    
    for s2k = 1:NumSequences
        S2K = ismember(1:NumSequences,1:s2k)';
        if Speller == 4
            S2K = repmat(S2K,2*NumTrials,1);
            trainD = TrainD;
            classlbl = ClassLabel(ChanLabel==1);
            trainD(ClassLabel(ChanLabel==1)==3,:) = [];
            classlbl(classlbl==3) = [];
        else
            S2K = repmat(S2K,NumStimCodes*NumTrials,1);
            trainD = TrainD;
            classlbl = ClassLabel(ChanLabel==1);
        end
        traind = trainD(S2K,:);
        rem_ind = sum(isnan(traind),2)>0;
        traind = traind(~rem_ind,:);
        classlbl = ClassLabel(ChanLabel==1);
        classlbl = classlbl(S2K);
        classlbl = classlbl(~rem_ind);
        if Cfier == 1
            %Stepwise Linear Regression
            INMODEL = zeros(1,size(AvgTrainD{s2k},2));
%             classlbl(classlbl==2) = -1; % To compare with Blankertz
%             (below)
            penter_v = .005; premove_v = .01;
            while (penter_v<.9 && sum(INMODEL)<8)
                penter_v = penter_v*2;
                premove_v = premove_v*2;
                if premove_v>1
                    penter_v=.95;
                    premove_v = .99;
                end
                if Avv == 0
                    [B, SE, PVAL, INMODEL] = stepwisefit(traind,...
                        classlbl,'penter',penter_v,'premove',premove_v,'display','off');
                
                    % SWLDA used by the Blankertz group -- same results.
%                     cd('BBCI SWLDA')
%                     %First row is assigned to -1 in the classifier, second
%                     %row to 1
%                     classlbl2(1,:) = classlbl == -1;
%                     classlbl2(2,:) = classlbl == 1;
%                     [C2,maxVar] = train_SWLDAmatlab(traind', classlbl2,...
%                         'PEntry',penter_v,'PRemoval',premove_v)
%                     cd('..')
                else
                    [B, SE, PVAL, INMODEL] = stepwisefit(AvgTrainD{s2k},...
                        AvgClassLabel(AvgChanLabel==1),'penter',penter_v,...
                        'premove',premove_v,'display','off');
                end
            end
            Classweight = -B(INMODEL==1);
            SigLocs = find(INMODEL==1);
            
            
            
            %         [~, keepind] = sort(PVAL);
            %         if length(keepind)>10
            %             keepind = keepind(1:10);
            %         end
            %         SigLocs = keepind;
            %         TrainD2 = TrainD(:,SigLocs);
            %         [B, SE, PVAL, INMODEL, STATS] = stepwisefit(TrainD2,ClassLabel(ChanLabel==1));
            %         Classweight = -B(INMODEL==1);
            %         SigLocs = SigLocs(INMODEL==1);
            
            
        elseif Cfier == 2
            %LDA
            %Classifier Built from individual trials (dont think this is what I
            %want to do)
            if Avv == 0
                try
                [class,err,post,logp,coeff] = classify(traind, traind, classlbl);
                catch
                    coeff(1,2).linear = [];
                end
            elseif Avv == 1
                [class,err,post,logp,coeff] = classify(AvgTrainD{s2k},AvgTrainD{s2k},AvgClassLabel(AvgChanLabel==1));
            end
            Classweight = coeff(1,2).linear;
            SigLocs = find(Classweight~=0);
            
            
        elseif Cfier == 3
            %         % Steve's LDA
            %         %This produces a crazy gamma matrix similar to the regular LDA, which
            %         %doesnt work.  Figure out how to make regular LDA output a useful gamma
            %         %matrix and then come back to this
            %
            %         if Avv == 0
            %             cd (['C:\Users\' pcusr '\Google Drive\Year 4\'...
            %                 'ALS Project\Part2 - P300\P300Classifier\LDA Code\Redo'])
            %             [Gam,z,SigLocs] = runLDA2([],TrainD,ClassLabel(ChanLabel==1),10);
            %             cd (['C:\Users\' pcusr '\Google Drive\Year 4\'...
            %                 'ALS Project\Part2 - P300\P300Classifier'])
            %         elseif Avv == 1
            %             cd (['C:\Users\' pcusr '\Google Drive\Year 4\'...
            %                 'ALS Project\Part2 - P300\P300Classifier\LDA Code\Redo'])
            %             [Gam,z,SigLocs] = runLDA2([],AvgTrainD,AvgClassLabel(AvgChanLabel==1),10);
            %             cd (['C:\Users\' pcusr '\Google Drive\Year 4\'...
            %                 'ALS Project\Part2 - P300\P300Classifier'])
            %         end
            %         %Sometimes the polarity is reversed on the output.  We always want
            %         %the first class label (p300) to be greater
            %         if mean(z(1:length(z)/2))-mean(z(length(z)/2+1:end))>0
            %             Classweight = Gam;
            %         else
            %             Classweight= -Gam;
            %         end
        elseif Cfier == 4 %LASSO
            if Avv == 0
                [B2, STATS2] =  lasso(traind,classlbl,'CV',5);
            elseif Avv == 1
                [B2, STATS2] =  lasso(AvgTrainD{s2k},AvgClassLabel(AvgChanLabel==1),'CV',5);
            end
            %Retain variables in the model corresponding to the sparsest
            %model within one standard error of the minimum cross validated
            %MSE
            Classweight = -B2(:,STATS2.Index1SE);
            SigLocs = find(Classweight~=0);
            Classweight = Classweight(SigLocs);
            
            lassoPlot(B2,STATS2,'PlotType','CV')
            
            
            %         Xplus = [ones(size(TrainD,1),1) TrainD];
            %         fitSparse = Xplus * [STATS2.Intercept(STATS2.Index1SE); B2(:,STATS2.Index1SE)];
            %         corr(fitSparse,ClassLabel(ChanLabel==1)-fitSparse)
            %         figure
            %         plot(fitSparse,ClassLabel(ChanLabel==1)-fitSparse,'o')
            
        elseif Cfier == 5 %Regularized LDA
            if Avv == 0
                Mdl = fitcdiscr(traind,classlbl,...
                    'SaveMemory','on');
                [err,gamma,delta,numpred] = cvshrink(Mdl,...
                    'NumGamma',25,'NumDelta',25,'Verbose',1);
                if Figures_On == 1
                    figure('Name','Error vs number of predictors for regularized LDA classifier');
                    plot(err,numpred,'k.')
                    xlabel('Error rate');
                    ylabel('Number of predictors');
                end
                low60 = min(min(err(numpred <= 60)));
                lownum = min(min(numpred(err == low60)));
                [p, q] = find((err == low60) & (numpred == lownum));
                p = p(1); q = q(1);  %% added this in the case of multiple minimums
                idx = sub2ind(size(delta),p,q);
                [gamma(p); delta(idx)]
                
                if Figures_On ==1
                    figure('Name','RegLDA')
                    plot(Mdl.Coeffs(1,2).Linear); hold on;
                    plot(Mdl.Coeffs(1,2).Linear,'r');
                end
                Mdl.Gamma = gamma(p);
                Mdl.Delta = delta(idx);
                Classweight = Mdl.Coeffs(1,2).Linear;
                SigLocs = find(Classweight~=0);
                Classweight = Classweight(SigLocs);
                
            elseif Avv == 1
                Mdl = fitcdiscr(AvgTrainD{s2k},AvgClassLabel(AvgChanLabel==1),...
                    'SaveMemory','on');
                [err,gamma,delta,numpred] = cvshrink(Mdl,...
                    'NumGamma',9,'NumDelta',9,'Verbose',1);
                if Figures_On==1
                    figure('Name','Error vs number of predictors for regularized LDA classifier');
                    plot(err,numpred,'k.')
                    xlabel('Error rate');
                    ylabel('Number of predictors');
                end
                low60 = min(min(err(numpred <= 60)));
                lownum = min(min(numpred(err == low60)));
                [p, q] = find((err == low60) & (numpred == lownum));
                p = p(1); q = q(1);  %% added this in the case of multiple minimums
                idx = sub2ind(size(delta),p,q);
                [gamma(p); delta(idx)]
                
                if Figures_On==1
                    figure('Name','RegLDA')
                    plot(Mdl.Coeffs(1,2).Linear); hold on;
                    plot(Mdl.Coeffs(1,2).Linear,'r');
                end
                Mdl.Gamma = gamma(p);
                Mdl.Delta = delta(idx);
                Classweight = Mdl.Coeffs(1,2).Linear;
                SigLocs = find(Classweight~=0);
                Classweight = Classweight(SigLocs);
            end
        elseif Cfier == 6
            %Manual classifier creation
            SigLocs = [];
            if s2k==5
                figure('Position',[100 100 1200,500])
                hold on;
                errorbar(nanmean(traind(classlbl==1,:)),nanstd(traind(classlbl==1,:))/sqrt(sum(classlbl==1)),'b')
                errorbar(nanmean(traind(classlbl==2,:)),nanstd(traind(classlbl==2,:))/sqrt(sum(classlbl==2)),'r')
                xlim([1 size(traind,2)]); title('Select 5 Points');
                
                [SigLocs_i,~] = ginput(20);
                SigLocs_i = round(SigLocs_i);
                Classweight = sign(nanmean(traind(classlbl==1,SigLocs_i))-nanmean(traind(classlbl==2,SigLocs_i)));
                
                
                [B, SE, PVAL, INMODEL] = stepwisefit(traind(:,SigLocs_i),...
                    classlbl,'display','off');
                Classweight = -B(INMODEL==1);
                SigLocs = SigLocs_i(find(INMODEL==1));
                
            end
            
        end
        
        %Construct Classifier Matrix
        CM{s2k} = zeros(length(SigLocs),4);
        if isempty(CM{s2k})
            disp('No significant features!');
            CM{s2k} = [];
            skiptest = 1;
        else
            CM{s2k}(:,1) = inpch(SigLocs); %Input Channel
            Ctms = mod(SigLocs,length(LDAtimes));  Ctms(Ctms==0)=length(LDAtimes);
            CM{s2k}(:,2) = LDAtimes(Ctms); %Input element (sample #)
            CM{s2k}(:,3) = 1; %Output channel
            CM{s2k}(:,4) = Classweight; %Weight
            % Changed from CM(:,4) = Classweight(SigLocs);  to reflect changes
            % to SW classifier that limits features to a specified number.
            
            
            
            
        end
        waitbar(.6+.3*s2k/NumSequences,wtbrr);
    end
    temp = cellfun(@(x) ~isempty(x),CM);
    bigCM = max(find(temp));
    if isempty(bigCM)
        cm_full = [];
    else
        cm_full = CM{bigCM};
    end
    
    
else %Testrun - do not build classifier, prompt for a previously built one
    skiptest = 0;
    %Load Classifier Data
    try
        disp(['Loading classifier: ' NameBase ClassTypes{Cfier} 'Classifier.txt']);
        CM = dlmread([NameBase ClassTypes{Cfier} 'Classifier.txt'],'\t');
    catch
        disp(['Loading classifier: ' NameBase2 ClassTypes{Cfier} 'Classifier.txt']);
        CM = dlmread([NameBase2 ClassTypes{Cfier} 'Classifier.txt'],'\t');
    end
    cm_full = CM;
    %Unnecessary for P300
    %     NormalizerGain = CM(end-1,1);
    %     NormalizerOffset = CM(end,1);
    %     CM = CM(1:end-2,:);
    
    % Name = 'AndrewGrid-gUSB';
    % CSession = '002';
    % CRun = '02';
    % CM = dlmread(['C:\Documents and Settings\amg5106\My Documents\'...
    %         'Year2\Research\BCI2000\BCI2000sourcetree\data\AndrewGrid-gUSB002\'...
    %         'AndrewGrid-gUSBS002R01_SWClassifier.txt'],'\t');
    
end


if FreeSpell == 0
    %Looking at the expected classification for trial averages
    for ii = 1:size(AvgTrainD{end},1)
        cursA(ii) = 0;
        for jj = 1:size(cm_full,1)
            ffr = (cm_full(jj,1)-1)*length(LDAtimes)+1:cm_full(jj,1)*length(LDAtimes);
            ffrx = LDAtimes==cm_full(jj,2);
            cursA(ii) = cursA(ii) + cm_full(jj,4)*AvgTrainD{end}(ii,ffr(ffrx));
            AvgTrainD{end}(ii,ffr(ffrx));
        end
    end
    %What is optimal decision boundary?
    pcurs = cursA(AvgClassLabel(AvgChanLabel==1)==1);
    ncurs = cursA(AvgClassLabel(AvgChanLabel==1)==2);
    k=1;
    rangeC = linspace(min(ncurs),max(pcurs),100);
    for thr = rangeC %Cycle through possible thresholds
        pp = sum(pcurs>thr)/length(pcurs);
        pn = sum(ncurs<thr)/length(ncurs);
        psums(k) = pp+pn; k=k+1;
    end
    [ma mi] = max(psums);
    %Set the gain and the offset of the Normalizer
    NormalizerOffset = rangeC(mi);
    NormalizerGain = 1/var(cursA);
    
    
    %Looking at the expected classification for individual segments
    for ii = 1:size(TrainD,1)
        curs(ii) = 0;
        for jj = 1:size(cm_full,1)
            ffr = (cm_full(jj,1)-1)*length(LDAtimes)+1:cm_full(jj,1)*length(LDAtimes);
            ffrx = LDAtimes==cm_full(jj,2);
            curs(ii) = curs(ii) + cm_full(jj,4)*TrainD(ii,ffr(ffrx));
        end
    end
    for tt = 1:NumTrials*NumStimCodes
        trlindex = (tt-1)*NumSequences+1:(tt)*NumSequences;
        cursA2(tt) = mean(curs(trlindex));
    end
    
    if Figures_On == 1
        figure('Name','Histogram of target and non-target amplitudes')
        subplot(211);
        hist(cursA(AvgClassLabel(AvgChanLabel==1)==1),10);
        hold on;
        hist(cursA(AvgClassLabel(AvgChanLabel==1)==2),10);
        h = findobj(gca,'Type','patch');
        set(h(1),'FaceColor','r','EdgeColor','w','FaceAlpha',.5);
        %weird that this labeling of h(1) and h(2) is backwards
        set(h(2),'FaceColor','b','EdgeColor','w','FaceAlpha',.5);
        line([NormalizerOffset NormalizerOffset],[0 20],'Color','k','LineWidth',3)
        subplot(212);
        hist(curs(ClassLabel(ChanLabel==1)==1),100);
        hold on;
        hist(curs(ClassLabel(ChanLabel==1)==2),100);
        h = findobj(gca,'Type','patch');
        set(h(1),'FaceColor','r','EdgeColor','w','FaceAlpha',.5);
        set(h(2),'FaceColor','b','EdgeColor','w','FaceAlpha',.5);
    end
    
    
    
    %Instead of the line on the histogram plot, do an ROC curve
    %If we consider "truth" as being a left trial
    k=1;
    for thr = rangeC
        tp(k) = sum(pcurs>thr); fn(k) = sum(pcurs<thr);
        tn(k) = sum(ncurs<thr); fp(k) = sum(ncurs>thr); k = k+1;
    end
    tpr = tp./(tp+fn);
    fpr = fp./(fp+tn);
    auckinda = tpr.*(1-fpr);
    [~, emaxi] = max(auckinda);
    rangeC(emaxi);
    lineslope = 1;
    rangeL = 0:.01:1;
    deg45line = lineslope*(rangeL-fpr(emaxi))+tpr(emaxi);
    if Figures_On == 1
        figure('Name','FPR vs TPR')
        plot(fpr,tpr); hold on
        plot(rangeL,deg45line,'r'); xlim([0 1]); ylim([0 1]);
    end
    %Do the two methods produce the same result?
    if mi==emaxi
        disp('Yay!');
    else
        disp('Different :(');
    end
end

%     CM(size(CM,1)+1,:) = [NormalizerGain 0 0 0];
%     CM(size(CM,1)+1,:) = [NormalizerOffset 0 0 0];

%% View Results
clear MinEvidence
for s2k = 1:NumSequences
    
    if BuildClassifier == 1
        cm = CM{s2k};
    else
        cm = CM;
    end
    
    if isempty(cm)
        
        WORD_full(s2k,:) = repmat({'/'},1,length(TTSpell));
        if FreeSpell == 0
            WORDD_full(s2k) = sum(cell2mat(WORD_full(s2k,:)) == TTSpell)/length(TTSpell)*100;
        end
        allsums = [];
        
    else
        S2K = ismember(1:NumSequences,1:s2k)';
        S2K = repmat(S2K,NumStimCodes*NumTrials,1);
        if FreeSpell == 0 %Find the classified choice and not choice data
            %Absolutely no difference in whether the data is classified using the
            %Averaged or non averaged data. The result of choice and notchoice are
            %the same.  FOr this reason we remove the functionality of switching
            %between Avv=0 and Avv=1
            
            
            
            
            %     if Avv == 0
            %         ChoiceOut = []; NotChoiceOut = [];
            %         % ChoiceOut = coeff(1,2).const*ones(sum(ClassLabel(ChanLabel==1)=='C'),1);
            %         % NotChoiceOut = coeff(1,2).const*ones(sum(ClassLabel(ChanLabel==1)=='N'),1);
            %         ChoiceL = ClassLabel(1:NumStimCodes*NumSequences*NumTrials)=='C';
            %         NChoiceL = ClassLabel(1:NumStimCodes*NumSequences*NumTrials)=='N';
            %         %Plot the classification results with the generated classifier
            %         for i = 1:size(ClassifierMatrix,1)
            %             tform = EEGData((ClassifierMatrix(i,1)-1)*(NumStimCodes*NumSequences*NumTrials)+1:...
            %                 ClassifierMatrix(i,1)*(NumStimCodes*NumSequences*NumTrials),ClassifierMatrix(i,2))...
            %                 *ClassifierMatrix(i,4);
            %             ChoiceOut = [ChoiceOut tform(ChoiceL)];
            %             NotChoiceOut = [NotChoiceOut tform(NChoiceL)];
            %         end
            %         %Consider for each trial what the sum of the choice and not choice
            %         %sequences are
            %         for i = 1: NumTrials
            %             choicesums{i} = sum(reshape(sum(ChoiceOut(i*size(targetstim,2)*NumSequences...
            %             -(size(targetstim,2)*NumSequences-1):i*size(targetstim,2)*NumSequences,:),2),...
            %             NumSequences,size(targetstim,2)),1);
            %
            %             notchoicesums{i} = sum(reshape(sum(NotChoiceOut(i*...
            %                 (NumStimCodes-size(targetstim,2))*NumSequences...
            %             -((NumStimCodes-size(targetstim,2))*NumSequences-1):i*...
            %             (NumStimCodes-size(targetstim,2))*NumSequences,:),2),...
            %             NumSequences,NumStimCodes-size(targetstim,2)),1);
            %
            %             figure
            %             hist(notchoicesums{i});hold on;
            %             hist(choicesums{i});
            %             h = findobj(gca,'Type','patch');
            %             set(h(1),'FaceColor','r');
            %         end
            %
            %         figure
            %         plot(range/fs,EEGData(ClassLabel=='N'&ChanLabel==Clbl,:)','r');
            %         hold on;
            %         plot(range/fs,EEGData(ClassLabel=='C'&ChanLabel==Clbl,:)','b');
            %
            %         figure
            %         errorbar(mean(EEGData(ClassLabel=='C'&ChanLabel==Clbl,range),1),...
            %             std(EEGData(ClassLabel=='C'&ChanLabel==Clbl,range),[],1)/...
            %             sqrt(size(EEGData(ClassLabel=='C'&ChanLabel==Clbl,range),1))); hold on;
            %         errorbar(mean(EEGData(ClassLabel=='N'&ChanLabel==Clbl,range),1),...
            %             std(EEGData(ClassLabel=='N'&ChanLabel==Clbl,range),[],1)/...
            %             sqrt(size(EEGData(ClassLabel=='N'&ChanLabel==Clbl,range),1)),'r');
            %         xlim([0 length(range)])
            %    else
            AllOut = [];
            AvgAllOut = [];
            
            
            
            % THIS IS DONE ON UNPROCESSED DATA -- NO ARTIFACTING IN
            % EEGDATA.  IF USING DATA WITH SOME TRIALS REMOVED, THIS METHOD
            % WOULD NOT WORK.  THIS IS NOT THAT PROBLEMATIC.  IT SHOWS WHAT
            % THE ACCURACY MIGHT BE LIKE ON REAL DATA WHERE REAL TIME
            % ARTIFACT REJECTION IS NOT OCCURRING.  THE CLASSIFIER IS STILL
            % BUILT WITH ARTIFACTED DATA.  WAIT NO THIS IS DONE ASSUMING
            % ALL FILES HAVE THE SAME NUMBER OF SEQUENCES
            
            %Plot the classification results with the generated classifier
            for i = 1:size(cm,1)
                tform = EEGData((cm(i,1)-1)*(NumStimCodes*NumTrials*NumSequences)+1:...
                    cm(i,1)*(NumStimCodes*NumTrials*NumSequences),cm(i,2)-range(1)+1)*cm(i,4);
                AllOut = [AllOut tform];
                avgtform = AvgEEGData{s2k}((cm(i,1)-1)*(NumStimCodes*NumTrials)+1:...
                    cm(i,1)*(NumStimCodes*NumTrials),cm(i,2)-range(1)+1)*cm(i,4);
                AvgAllOut = [AvgAllOut avgtform];
            end
            clear allsums choicesums notchoicesums avg_choicesums ...
                avg_notchoicesums
            for i = 1: NumTrials
                %take sequences from a single trial
                trl_ind = (i-1)*NumStimCodes*NumSequences+1:i*NumStimCodes*NumSequences;
                avg_trl_ind = (i-1)*NumStimCodes+1:i*NumStimCodes;
                %only keep a certain number of sequences
                tmp_ind = trl_ind(S2K(trl_ind));
                %sum the data features
                temp_sum = nansum(AllOut(tmp_ind,:),2)';
                avg_allsums{i} = nansum(AvgAllOut(avg_trl_ind,:),2)';
                %reshape by number of sequences into stim codes
                temp_sum = reshape(temp_sum,s2k,NumStimCodes);
                if size(temp_sum,1)>1
                    allsums{i} = sum(temp_sum);
                else
                    allsums{i} = temp_sum;
                end
                
                cclabel = find(ismember(1:NumStimCodes,targetstim(i,:))==1);
                nnlabel = find(~ismember(1:NumStimCodes,targetstim(i,:))==1);
                choicesums{i} = allsums{i}(1,cclabel)/s2k;
                notchoicesums{i} = allsums{i}(1,nnlabel)/s2k;
                avg_choicesums{i} = avg_allsums{i}(1,cclabel);
                avg_notchoicesums{i} = avg_allsums{i}(1,nnlabel);
                
                
                MinEvTarget = cellfun(@(x) sum(x),avg_choicesums);
                MinEvNonTarget = cellfun(@(x,y) max(x)+max(y),avg_choicesums,avg_notchoicesums);
                MinEvidence(s2k) = mean(MinEvTarget-MinEvNonTarget);
                
                
                %         figure
                %         hist(notchoicesums{i});hold on;
                %         hist(choicesums{i});
                %         h = findobj(gca,'Type','patch');
                %         set(h(1),'FaceColor','r');
                %         if Figures_On == 1
                %             warning('off','all')
                %             figure
                %             scatter(cclabel,choicesums{i},'r');
                %             hold on
                %             scatter(nnlabel,notchoicesums{i},'b');
                %             for scc = 1:NumStimCodes
                %                 if Speller==1
                %                     inR = strcmp(num2str(scc),SR{1,i});
                %                     if sum(inR)>0
                %                         incodesR = find(strcmp(num2str(scc),SR{1,i}));
                %                         for icR = 1:length(incodesR)
                %                             text(double(scc),icR*2,TargSymb{incodesR(icR)})
                %                         end
                %                     else
                %                         incodesC = find(strcmp(num2str(scc),SC{1,i}));
                %                         for icC = 1:length(incodesC)
                %                             text(double(scc),icC*2,TargSymb{incodesC(icC)})
                %                         end
                %                     end
                %                 elseif Speller == 2 || Speller == 3
                %                     inR = scc==SR{i};
                %                     if sum(inR)>0
                %                         incodesR = find(scc==SR{i});
                %                         for icR = 1:length(incodesR)
                %                             text(double(scc),icR*2,TargSymb{incodesR(icR)})
                %                         end
                %                     else
                %                         incodesC = find(scc==SC{i});
                %                         for icC = 1:length(incodesC)
                %                             text(double(scc),icC*2,TargSymb{incodesC(icC)})
                %                         end
                %                     end
                %                 end
                %             end
                %             ylim([min([choicesums{i} notchoicesums{i} 2]) ...
                %                 max([choicesums{i} notchoicesums{i} ...
                %                 2*max([length(incodesR) length(incodesC)])])]);
                %             warning('on','all')
                %         end
            end
            
            
            %     if Figures_On == 1
            %         figure('Name','Averages of each stimulus code')
            %         for Clbl = 1:NumChans
            %             subplot(4,5,Clbl);
            %             plot(range/fs,AvgEEGData(AvgClassLabel==2&AvgChanLabel==Clbl,:)','r');
            %             hold on;
            %             plot(range/fs,AvgEEGData(AvgClassLabel==1&AvgChanLabel==Clbl,:)','b');
            %         end
            %
            %         figure('Name','Grand Averages with error bars')
            %         for Clbl = 1:NumChans
            %             subplot(4,5,Clbl); errorbar(range/fs,mean(AvgEEGData(AvgClassLabel==1&AvgChanLabel==Clbl,:),1),...
            %                 std(AvgEEGData(AvgClassLabel==1&AvgChanLabel==Clbl,:),[],1)/...
            %                 sqrt(size(AvgEEGData(AvgClassLabel==1&AvgChanLabel==Clbl,:),1))); hold on;
            %             errorbar(range/fs,mean(AvgEEGData(AvgClassLabel==2&AvgChanLabel==Clbl,:),1),...
            %                 std(AvgEEGData(AvgClassLabel==2&AvgChanLabel==Clbl,:),[],1)/...
            %                 sqrt(size(AvgEEGData(AvgClassLabel==2&AvgChanLabel==Clbl,:),1)),'r');
            %         end
            %     end
            %     figure
            %     hist(sum(NotChoiceOut,2),20);hold on;
            %     hist(sum(ChoiceOut,2),20);
            %     h = findobj(gca,'Type','patch');
            %     set(h(1),'FaceColor','r');
            %
        else % In Free Spelling mode
            
            %     if Avv == 0
            %         ChoiceOut = [];
            %         for i = 1:size(ClassifierMatrix,1)
            %             tform = EEGData((ClassifierMatrix(i,1)-1)*(NumStimCodes*NumSequences*NumTrials)+1:...
            %                 ClassifierMatrix(i,1)*(NumStimCodes*NumSequences*NumTrials),ClassifierMatrix(i,2))...
            %                 *ClassifierMatrix(i,4);
            %             ChoiceOut = [ChoiceOut tform];
            %         end
            %         for i = 1: NumTrials
            %             choicesums{i} = sum(reshape(sum(ChoiceOut(i*NumStimCodes*NumSequences...
            %                 -(NumStimCodes*NumSequences-1):i*NumStimCodes*NumSequences,:),2),...
            %                 NumSequences,NumStimCodes),1);
            %             notchoicesums{i} = [];
            %         end
            %     else
            
            %%%%%%%This needs to be updated to use average data.
            %             AllOut = [];
            %             for i = 1:size(cm,1)
            %
            %                 tform = EEGData((cm(i,1)-1)*(NumStimCodes*NumTrials*NumSequences)+1:...
            %                     cm(i,1)*(NumStimCodes*NumTrials*NumSequences),cm(i,2)-range(1)+1)*cm(i,4);
            %
            %                 AllOut = [AllOut tform];
            %
            %             end
            %
            %             clear choicesums notchoicesums
            %             for i = 1: NumTrials
            %                 %take sequences from a single trial
            %                 trl_ind = (i-1)*NumStimCodes*NumSequences+1:i*NumStimCodes*NumSequences;
            %                 %only keep a certain number of sequences
            %                 tmp_ind = trl_ind(S2K(trl_ind));
            %                 %sum the data features
            %                 temp_sum = sum(AllOut(tmp_ind,:),2)';
            %                 %reshape by number of sequences into stim codes
            %                 temp_sum = reshape(temp_sum,s2k,NumStimCodes);
            %                 if size(temp_sum,1)>1
            %                     choicesums{i} = sum(temp_sum);
            %                 else
            %                     choicesums{i} = temp_sum;
            %                 end
            %                 notchoicesums{i} = [];
            %             end
        end
        % end
        
        
        
        %% Trial Choice
        for i = 1:NumTrials
            %Sort the individual stimulus code sums
            
            [TrialVal(i,:) TrialChoice(i,:)] = sort([avg_choicesums{i} avg_notchoicesums{i}],'descend');
            
            %Find the actual stimulus codes associated with these
            trlcode = [targetstim(i,:) find(~ismember(1:NumStimCodes,targetstim(i,:))==1)];
            
            ll=1;
            clear trlchc chcsums
            for jj = 1:size(TrialChoice,2) %Loop through code 1
                for kk = (jj+1):size(TrialChoice,2) %Loop through code 2
                    %For each pair of codes, find the pairing of actual codes
                    trlchc(ll,:) = trlcode(TrialChoice(i,[jj kk]));
                    %Find the sum of these codes
                    chcsums(ll,:) = sum(TrialVal(i,[jj,kk]));
                    ll = ll+1;
                end
            end
            [~, chcsumsort] = sort(chcsums,1,'descend');
            
            matchfind = 0;
            while matchfind == 0
                for jj = chcsumsort'
                    jj;
                    
                    if Speller == 1
                        if ismember(trlchc(jj,1),cellfun(@str2num,SR{i}))
                            mtch = find(cellfun(@str2num,SR{i})==trlchc(jj,1)&...
                                cellfun(@str2num,SC{i})==trlchc(jj,2));
                        else
                            mtch = find(cellfun(@str2num,SC{i})==trlchc(jj,1)&...
                                cellfun(@str2num,SR{i})==trlchc(jj,2));
                        end
                    elseif Speller == 2 || Speller == 3  %DOUBLE CHECK THIS FOR RSVP CONDITION--- SHOULD BE IN SP4 CONDITION, NO?
                        if ismember(trlchc(jj,1),SR{i})
                            mtch = find(SR{i}==trlchc(jj,1)&...
                                SC{i}==trlchc(jj,2));
                        else
                            mtch = find(SC{i}==trlchc(jj,1)&...
                                SR{i}==trlchc(jj,2));
                        end
                    elseif Speller == 4
                        %Separated SP4 from SP2 and 3 --- 6/8/17
                        mtch = trlcode(TrialChoice(i,1));
                    end
                    
                    if ~isempty(mtch)
                        matchfind = 1;
                        StimChoice(i,:) = trlchc(jj,:);
                        TargChoice(i) = mtch;
                        
                        break;
                    end
                end
            end
            
            
        end
        
        idx = sub2ind(size(TargSymb),AllTrials,TargChoice');
        WORD_full(s2k,:) = TargSymb(idx);
        if FreeSpell == 0
            try
                corrword = 0;
                for itmp = 1:size(WORD_full,2)
                    if strcmp(WORD_full{s2k,itmp},GoalLetterF{itmp});
                        corrword = corrword+1;
                    end
                end
                WORDD_full(s2k) = corrword/length(GoalTrialF)*100
            catch
                disp('not calculating accuracy correctly')
                WORDD_full(s2k) = NaN;
            end
        end
    end
end


waitbar(.9,wtbrr);

SequencesUsed = (1:NumSequences)';
Guess = WORD_full;
Accuracy = WORDD_full';
if BuildClassifier == 1
    NumFeatsUsed = cellfun(@(x) size(x,1),CM)'
else
    NumFeatsUsed = repmat([num2str(size(CM,1)) '*'],NumSequences,1);
end
if length(SequencesUsed)==1
    SequencesUsed
    Guess
    Accuracy
    NumFeatsUsed
else
    table(SequencesUsed, Guess, Accuracy, NumFeatsUsed)
end
Accuracy'
WORD_Acc = Accuracy(end);
%% Save parameter files
if BuildClassifier == 1
    try
        cffile = [bcifold '\data\' Name Sessions{max_sess_i} '\'...
            Name 'S' Sessions{max_sess_i} 'R' NewRunName '_' ClassTypes{Cfier} 'Classifier.txt'];
        
    catch
        cffile = ['P:\ALS Proj Data\' Name Sessions{max_sess_i} '\'...
            Name 'S' Sessions{max_sess_i} 'R' NewRunName '_' ClassTypes{Cfier} 'Classifier.txt'];
    end
    dlmwrite(cffile,cm_full,'delimiter','\t');
    skiptest = 0;
    
    
    %save in parameter file
    prmcls = num2str(reshape(cm_full',1,size(cm_full,1)*size(cm_full,2)));
    %Change Subject Name
    fpm2{1}{SNInd} = [fpm2{1}{SNInd}(1:SBegInd+1) Name...
        fpm2{1}{SNInd}(SEndInd-1:end)];
    %Change Session
    fpm2{1}{SEInd} = [fpm2{1}{SEInd}(1:SEBegInd+1) Sessions{max_sess_i}...
        fpm2{1}{SEInd}(SEEndInd-1:end)];
    %Change Spatial Filter
    prmsf = [num2str(size(SpatFilt))  ' ' num2str(reshape(SpatFilt',1,...
        size(SpatFilt,1)*size(SpatFilt,2)))];
    fpm2{1}{SFInd} = [fpm2{1}{SFInd}(1:SFBegInd+1) prmsf...
        fpm2{1}{SFInd}(SFEndInd-1:end)];
    %Change Classifier
    fpm2{1}{LinInd} = [fpm2{1}{LinInd}(1:BegInd2+1) num2str(size(cm_full,1)) ...
        fpm2{1}{LinInd}(EndInd2:BegInd+1) prmcls...
        fpm2{1}{LinInd}(EndInd-1:end)];
    %Change Sensitivity
    if ~isempty(STVYInd) & isfield(c,'Sensitivity')
        fpm2{1}{STVYInd} =  [fpm2{1}{STVYInd}(1:STVYBegInd+1) num2str(c.Sensitivity.NumericValue)  ' 1 1 6 ' ...
            fpm2{1}{STVYInd}(STVYEndInd:end)];
    end
    
    
    %WordBank
    
    %Change the Word to be spelled
    if ctype == 1
        WordBank = {'JAWS','LINK','PARK','TEN8','C3PO','R2D2','CART','MIND',...
            'ALOE','PINK','BLUE','GREY','TAPE','CAPE','BURN','WIND','WILD',...
            'JUMP','PUMP','CAMP','MAPS','BARN','COAL','BOOT','TOLL','HUNT',...
            'FREE','COOL','MARK','DING','TINY','ZOOM','ZONE','XRAY','EXIT',...
            'NOUN','BOUT','BLIP','WALK','WORM','WISE'};
        temp = randperm(length(WordBank));
        fpm2{1}{TTSInd} = [fpm2{1}{TTSInd}(1:TTSBegInd+1) WordBank{temp(1)} ...
            fpm2{1}{TTSInd}(TTSEndInd-1:end)];
    elseif ctype == 2
        temp = randperm(length(SymbSymb(1:end-5)));
        temp = cat(2,SymbSymb{temp(1:4)});
        fpm2{1}{TTSInd} = [fpm2{1}{TTSInd}(1:TTSBegInd+1) temp...
            fpm2{1}{TTSInd}(TTSEndInd-1:end)];
    elseif ctype == 5
        temp = randperm(4);
        temp2 = [1 1 2 2];
        fpm2{1}{TTSInd} = [fpm2{1}{TTSInd}(1:TTSBegInd+1) [num2str([4 temp2(temp) 1 1]) ' %']...
            fpm2{1}{TTSInd}(TTSEndInd-1:end)];
    else
        fpm2{1}{TTSInd} = [fpm2{1}{TTSInd}(1:TTSBegInd+1) '% ' ...
            fpm2{1}{TTSInd}(TTSEndInd-1:end)];
    end
    %Change the letters to be displayed
    if ctype ~= 4
        fpm2{1}{DispResInd} = [fpm2{1}{DispResInd}(1:DispResBegInd+1) '1 1 0 1 ' ...
            fpm2{1}{DispResInd}(DispResEndInd:end)];
    end
    %         %Add Minimum Evidence --- NOT DOING THIS YET
    %         if ~isempty(AccEvInd)
    %             if MinEvidence(bigCM) < 0
    %                 disp('Target evidence not strong enough to set MinEvidence parameter')
    %                 MinEvidence(bigCM) = 0;
    %             end
    %             prmme = [num2str(MinEvidence(bigCM)) ' 0 0 % '];
    %             fpm2{1}{MinEvInd} = [fpm2{1}{MinEvInd}(1:MBegInd+1) prmme...
    %                 fpm2{1}{MinEvInd}(MEndInd-1:end)];
    %         end
    if ~isempty(ImpInd)
        %Change Acquisition Mode from Impedance to Analog Acquisition
        fpm2{1}{ImpInd} = [fpm2{1}{ImpInd}(1:ImpBegInd+1) '0 0 0 2 ' ...
            fpm2{1}{ImpInd}(ImpEndInd:end)];
    end
    
    
    fpm2f = fpm2;  %Parameter for FreeSpelling
    temp = regexp(fpm2f{1}{IntInd},'=');
    fpm2f{1}{IntInd} = [fpm2f{1}{IntInd}(1:temp) ' 1 ' fpm2f{1}{IntInd}(temp+4:end)]; %FreeSpell
    fpm2f{1}{TTSInd} = []; %No text to spell
    
    %Save new parameter file
    %Here, should save parameter files for use in the trainer,
    %symbolic, dual, and notepad spellers
    prmfile = [bcifold '\data\' Name Sessions{max_sess_i} '\' ptype '_'...
        Name 'S' Sessions{max_sess_i} 'R' NewRunName '.prm'];
    prmfileshort = ['...\\data\\' Name Sessions{max_sess_i} '\\' ptype '\_'...
        Name 'S' Sessions{max_sess_i} 'R' NewRunName '.prm'];
    
    if dmmy == 1 %save to dummy file
        dlmcell([bcifold '\data\' Name Sessions{max_sess_i} '\dummy_'...
            Name 'S' Sessions{max_sess_i} 'R' NewRunName '.prm'],fpm2{1});
    else
        dlmcell([bcifold '\data\' Name Sessions{max_sess_i} '\' ptype '_'...
            Name 'S' Sessions{max_sess_i} 'R' NewRunName '.prm'],fpm2{1});
    end
    
    
end

close(wtbrr)
end

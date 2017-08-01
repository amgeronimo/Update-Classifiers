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


function [CM, acc_fromCfier, traind, classlbl, ChannelNames, freqs, prmfile, eyye] = ...
    BuildMIClassifier_teleBCI2(Name,Sessions,Runs,BuildClassifier,Cfier,...
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
Cursor = [];
ResCode = [];
TargCode = [];
Feedback = [];
NewRunName = [];
GazeX = [];
GazeY = [];
SourceTime = [];
StimulusTime = [];
TTSpell = [];
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
                ChannelNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};
                EEGloc = 1:8; EOGloc =1;
            else
                capt = 'siggen';
                ChannelNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};
                EEGloc = 1:8; EOGloc = 1; HBF = 0;
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
        
        
        
        
        
        
        
        fs = c.SamplingRate.NumericValue;
        sbs = c.SampleBlockSize.NumericValue;
        %Trim Data and StimulusCodes
        trimMax = round((c.PreRunDuration.NumericValue-.5)*fs);
        Data = [Data;  a(trimMax:end,:)];
        Cursor = [Cursor; double(b.CursorPosX(trimMax:end))];
        TargCode = [TargCode; double(b.TargetCode(trimMax:end))];  %?? tempTC =  tempTC(1:c.PreRunDuration.NumericValue*fs+9)=0;
        ResCode = [ResCode; double(b.ResultCode(trimMax:end))];
        Feedback = [Feedback; double(b.Feedback(trimMax:end))];
        SourceTime = [SourceTime; b.SourceTime(trimMax:end)];
        StimulusTime = [StimulusTime; b.StimulusTime(trimMax:end)];
        if isfield(b,'EyetrackerLeftEyeGazeX')
            GazeX = [GazeX; b.EyetrackerLeftEyeGazeX(trimMax:end)];
            GazeY = [GazeY; b.EyetrackerLeftEyeGazeY(trimMax:end)];
        else
            GazeX = [GazeX; NaN(size(a(trimMax:end,:),1),1)];
            GazeY = [GazeY; NaN(size(a(trimMax:end,:),1),1)];
        end
        
        NewRunName = strcat(NewRunName,Runs{ss}{rr})
        
        %Load Classifier
        ClassifierUsed{kk} = c.Classifier.NumericValue;
        ClassifierUsed{kk}(:,4) = cellfun(@str2num,c.Classifier.Value(:,4));
        
        %Number of channels
        NumChans{kk} = size(a,2);
        
        
        
        %         %Duration of the Stimulus State -- This is how long the stimulus is on
        %         StateDuration = c.StimulusDuration.NumericValue*fs/1000;
        
        %Number of trials
        NumTrials{kk} = c.NumberOfTrials.NumericValue;
        
        if isfield(b,'EyetrackerLeftEyeGazeX')
            AppPos = [AppPos; repmat([c.WindowLeft.NumericValue c.WindowTop.NumericValue ...
                c.WindowWidth.NumericValue c.WindowHeight.NumericValue],NumTrials{kk},1)];
        end
        
        %Assign each trial to a run number
        AllTrials = [AllTrials; kk*ones(NumTrials{kk},1)];
        
        %Load Spatial Filter
        SpatFiltUsed{kk} = c.SpatialFilter.NumericValue;
        
        %     sdtemp = diff(double(TargCode));
        %     tstart = find(sdtemp>0);
        %     tend = find(sdtemp<0);
        %     if length(tstart) ~= sum(cell2mat(NumTrials))
        % disp('run did not complete')
        %     end
        
        kk = kk+1;
    end
end
Data = Data';

disp(['Time to load data ' num2str(toc)])
waitbar(.05,wtbrr);
%  range = (1:round(fs));

%% Check Variables
%Check if free spelling
% FreeSpell = sum(StimulusType)==0&sum(StimulusCode)~=0;

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
    ClassifierUsed = ClassifierUsed{1};
end


%If using the CNEamp - remove the first 2 channels
if isfield(c,'ComPort')==1
    Data = Data(3:end,:);
    cneamp=1;
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
    plot(TargCode(1:sbs:end),'g');
    plot(-BadTime*10,'k')
end
%figure
%plot(Data(:,1)*1000); hold on; plot(b.SourceTime,'r');



%% Artifact Correction
TargCodeUNART = TargCode;

%Ocular artifacts
if isempty(EOGloc)
    disp('Cannot reject artifacts because no EOG channels were recorded');
    Art = 0;
end
if Art == 1 %Remove data displaying artifacts
    disp('Removing dead data....');
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
        TargCode(logical(SmLoc))=0; %= StimulusCode(~ArtLoc);
    else
        SmLoc = [];
    end
    
    disp('Rmoving superlarge data....');
    
    %Remove superlarge data artifacts 1s before, 5 seconds after
    LInt = find(abs(Data(EOGloc(1),:))>1e3);
    if ~isempty(LInt)
        Lrange = fix(-1*fs:5*fs);
        LInt = repmat(LInt,length(Lrange),1)+repmat(Lrange,length(LInt),1)';
        LInt = LInt(:);
        LInt = LInt(LInt>0 & LInt<size(Data,2));
        LLoc = zeros(size(Data,2),1);
        LLoc(LInt)=1;
        TargCode(logical(LLoc))=0; %= StimulusCode(~ArtLoc);
    else
        LLoc = [];
    end
    
     disp('Removing ocular artifacts....');
    
    ArtInt = find(Data(EOGloc(1),:)>50&circshift(Data(EOGloc(1),:),[0 1])<50);
    ArtInt = [ArtInt find(Data(EOGloc(1),:)<-50&circshift(Data(EOGloc(1),:),[0 1])>50)];
    
    Arange = fix(-.1*fs:1*fs);
    ArtInt = repmat(ArtInt,length(Arange),1)+repmat(Arange,length(ArtInt),1)';
    %     ArtInt = ArtInt';
    ArtInt = ArtInt(:);
    ArtInt = ArtInt(ArtInt>0 & ArtInt<size(Data,2));
    ArtLoc = zeros(size(Data,2),1);
    ArtLoc(ArtInt)=1;
    
    
    TargCode(logical(ArtLoc))=0; %= StimulusCode(~ArtLoc);
    
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
        plot(TargCode,'k');
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
    TargCodeUNART = TargCode;
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
%by filtering out low frequency data and performing peak detection
if HBF == 1
    [bb,aa] = butter(4,[3 60]/(fs/2));
    HBdata = filtfilt(bb,aa,Data')';
    
    % figure
    % kk = 1;
    % for i = 3:5
    %     xx(kk) = subplot(3,1,kk); plot(Data(i,:)); hold on; plot(HBdata(i,:));
    %     kk=kk+1;
    % end
    % linkaxes(xx)
    
    temp = mean(HBdata(6:8,:));
    [maxt, mint] = peakdet(temp,5*std(temp));
    if ~isempty(maxt)
        maxt(maxt(:,1)<25|maxt(:,1)>size(Data,2)-25,:) = [];
        bpm = fs*60./diff(maxt(:,1));
        %remove heartbeats greater than 100 bpm
        keepbeep = bpm<=100;
        bpm(~keepbeep,:) = [];
        statschan = [mean(bpm) std(bpm)];
        if Figures_On==1
            figure
            subplot(1,4,1:3); plot(temp); hold on;
            scatter(maxt(keepbeep,1),maxt(keepbeep,2));
            subplot(2,4,4); hist(bpm,200);
        end
    else
        statschan = [0 100];
    end
    
    if statschan(2)<15 && statschan(1)>40
        disp(['Strong heartbeat at ' num2str(statschan(1)) ' bpm, proceeding removal of heartbeat'])
        
        
        Data2 = Data;
        tempind = repmat(maxt(keepbeep,1),1,50)+repmat(-24:25,sum(keepbeep),1);
        wnd = window(@hanning,50);
        for i = 1:8
            tempd = reshape(HBdata(i,tempind'),50,sum(keepbeep));
            %remove artifactual epochs
            bade = abs(mean(tempd))>5*abs(mean(mean(tempd)));
            tempd(:,bade) = [];
            signature(:,i) = wnd.*mean(tempd,2);
            
            fkb = find(keepbeep);
            for j = 1:length(fkb)
                Data2(i,maxt(fkb(j),1)+[-24:25]) = Data2(i,maxt(fkb(j),1)+[-24:25])-signature(:,i)';
            end
        end
        
        if Figures_On==1
            subplot(2,4,8); plot(signature);
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
            xx(1) = subplot(211); plot(Data'); hold on; plot(TargCodeUNART,'k','LineWidth',1.5);
            xx(2) = subplot(212); plot(Data2'); hold on; plot(TargCode,'k','LineWidth',1.5);
            linkaxes(xx);
        end
        Data = Data2;
    else
        disp('Heartbeat not detected')
    end
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
        %         DataF = Data*SpatFilt;f
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


LC = TargCode==1;
RC = TargCode==2;
LCt = LC-circshift(LC,[1 0]);
RCt = RC-circshift(RC,[1 0]);
if RCt(1)==-1
    RCt = circshift(RCt,[-1 0]);
end
if LCt(1)==-1
    LCt = circshift(LCt,[-1 0]);
end
if Figures_On==1
    figure
    mm(1)=subplot(311);plot(LCt);
    mm(2)=subplot(312);plot(RCt,'r');
    linkaxes(mm)
end

L_trialstart = find(LCt==1);
L_trialend = find(LCt==-1);
R_trialstart = find(RCt==1);
R_trialend = find(RCt==-1);

if BuildClassifier == 1 %added 11/1/13 -- only need to make sure there is enough data when building classifier.
    incL = find(L_trialend-L_trialstart<768); %make sure there is more than three seconds of trial
    L_trialstart = L_trialstart(~ismember(1:length(L_trialstart),incL));
    L_trialend = L_trialend(~ismember(1:length(L_trialend),incL));
    
    incR = find(R_trialend-R_trialstart<768);
    R_trialstart = R_trialstart(~ismember(1:length(R_trialstart),incR));
    R_trialend = R_trialend(~ismember(1:length(R_trialend),incR));
end

NTl = length(L_trialstart);
NTr = length(R_trialstart);




disp(['Done with myfile ' num2str(toc)]);
waitbar(.1,wtbrr);


%% Do CSP
if DoCSP == 1
    if BuildClassifier == 1
        Data_C = (RefFilt*ArtWeights2*Data)';
        [W] = CSP_SMR(Data_C,StimulusCodeUNArt,StimulusType,...
            NumChans, NumTrials, NumStimCodes,range,StateDuration);
        disp(['Done with CSP ' num2str(toc)]);
    elseif BuildClassifier == 0
    end
elseif DoCSP == 0
    if BuildClassifier == 1
        W = eye(size(RefFilt,1));
    else
    end
    disp('CSP not done');
end


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
%Can include this later -- for now just ignore eye tracking for audio speller
%
if sum(~isnan(GazeX))~=0 && sum(GazeX)~=0
    EyeTracking_inMIClassifier(GazeX,GazeY, {[L_trialstart...
        L_trialend], [R_trialstart R_trialend]})
end
%         [gazeacc, gazevar, gazeinvalid] = P300_EyeTracking_inP300Classifier_v3(...
%             GazeX,GazeY, Data(EOGloc,:), SC,SR, StimulusCodeUNArt, ...
%             SequencePhase, StimulusTypeUNArt, StateDuration, ...
%             NumTrials,NumSequences, Speller,...
%             fs, c, Figures_On, 0, TTSpell, AppPos);
%         disp(['Finished with eye tracking ' num2str(toc)]);
%     else
%         disp('No eye tracking data detected');
%     end
%
%     waitbar(.25,wtbrr);
%     eyye.gazeacc = gazeacc;
%     eyye.gazevar = gazevar;
%     eyye.gazeinvalid = gazeinvalid;

eyye =[];


%% Calculate Spectra
%The results of PowSpect match the workings of FFTFilter.cpp in the BCI2000
%distribution.  in ::Process(), put code
%bciout << "Output(" << i << ",0) = " << Output(i,0) << endl;
%This allows us to see the result of the Spectra Calculation.  in BCI2000,
%the input to the fft builds up over time, while the spectrogram function
%matlab starts with a full window of data, and then shifts.  Therefore, the
%output of PowSpect(1,1,1) should match the (wndw/8)th iteration of
%Output(0,0).

NumTrials = sum(cell2mat(NumTrials));

wndw = c.WindowLength.NumericValue*fs;

hm = hamming(wndw*2+1);
hm = hm(wndw+1:end-1);
PowSpect = single(zeros(ceil(wndw/2+.5),floor((size(DD,2)-wndw)/sbs)+1,NumChans));
for ch = 1:NumChans
    ch
    %Spectrogram computes the short-time fourier transform with data
    %transformed by the half hamming
    [sS fF tT ~] = spectrogram(DD(ch,:),hm,wndw-sbs,wndw,fs);
    
    %fft data is reorganized like it is online
    if mod(wndw,2)==0
    sS2 = [real(sS); flipud(imag(sS(2:end-1,:)))];
    else
    sS2 = [real(sS); flipud(imag(sS(2:end,:)))];  
    end
    
    %Then the power spectrum is computed by squaring and adding real and
    %imaginary parts
    normfactor = 1/(wndw);
    mxIndx = ceil(wndw/2+.5);
    clear PS;
    PS(1,:) = sS2(1,:).^2*normfactor;
    for i = 1:mxIndx-1
        PS(i+1,:) = (sS2(i+1,:).^2+sS2(wndw+1-i,:).^2)*normfactor;
    end
    if mod(wndw,2)==0
        PS(mxIndx,:)= sS2(mxIndx,:).^2*normfactor;
    end
    PowSpect(:,:,ch) = PS;
end
clear sS sS2 PS

disp(['Time for calculating spectra ' num2str(toc)]);
waitbar(.3,wtbrr);



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





% %Plot averages before and after epoch removal
% tic
% if Figures_On == 1
%     if Avv == 1
%         if FreeSpell == 0
%             figure('Name','Target and non-target ERPs of filtered data, before and after epoch removal')
%             for ch = 1:NumChans
%                 subplot(4,5,ch);
%                 plot(range/fs,mean(AvgEEGDataRR{ns}(AvgClassLabel==2&AvgChanLabel==ch,:)),'Color',[1 .7 .7]);
%                 hold on;
%                 plot(range/fs,mean(AvgEEGDataRR{ns}(AvgClassLabel==1&AvgChanLabel==ch,:)),'Color',[.7 .7 1]);
%                 plot(range/fs,mean(AvgEEGData{ns}(AvgClassLabel==2&AvgChanLabel==ch,:)),'Color','r');
%                 plot(range/fs,mean(AvgEEGData{ns}(AvgClassLabel==1&AvgChanLabel==ch,:)),'Color','b');
%                 title(ChannelNames{ch});
%             end
%         elseif FreeSpell == 1
%             figure('Name','Averaged ERPs of filtered data')
%             for ch = 1:NumChans
%                 subplot(4,5,ch);
%                 plot(range,AvgEEGDataRR{ns}(AvgChanLabel==ch,:));
%             end
%         end
%     else %Avv=0
%         if FreeSpell == 0
%             figure('Name','Target and non-target ERPs of filtered data, before and after epoch removal')
%
%             t1=sum(isnan(EEGData),2)>0;
%             t2=sum(isnan(EEGDataRR),2)>0;
%             badremoved = zeros(size(EEGDataRR));
%             badremoved(t1,:) = EEGDataRR(t1,:);
%
%
%             for ch = 1:NumChans
%                 subplot(4,5,ch);
%                 plot(range/fs,nanmean(EEGDataRR(ClassLabel==2&ChanLabel==ch,:)),'Color',[1 .7 .7]);
%                 hold on;
%                 plot(range/fs,nanmean(EEGDataRR(ClassLabel==1&ChanLabel==ch,:)),'Color',[.7 .7 1]);
%                 plot(range/fs,nanmean(badremoved(ClassLabel==2&ChanLabel==ch,:)),'Color',[.7 1 .7]);
%                 plot(range/fs,nanmean(EEGData(ClassLabel==2&ChanLabel==ch,:)),'Color','r');
%                 plot(range/fs,nanmean(EEGData(ClassLabel==1&ChanLabel==ch,:)),'Color','b');
%                 plot(range/fs,nanmean(badremoved(ClassLabel==1&ChanLabel==ch,:)),'Color','g');
%
%                 title(ChannelNames{ch});
%             end
%         elseif FreeSpell == 1
%             figure('Name','Averaged ERPs of filtered data')
%             for ch = 1:NumChans
%                 subplot(4,5,ch);
%                 plot(range,nanmean(EEGDataRR(ChanLabel==ch,:)),'b');
%             end
%         end
%     end
% end
%
% clear EEGDataRR AvgEEGDataRR
%
%
% disp(['Plot after epoch removal ' num2str(toc)]);


%% Sanity Check (3)
if Figures_On == 1
    rtms = find(RC)/fs;
    rtms = ismember(round(tT,2),round(rtms,2));
    ltms = find(LC)/fs;
    ltms = ismember(round(tT,2),round(ltms,2));
    figure('Name','Average target and non-target ERPs of filtered, down sampled, and epoch rejected data');
    for ch = 1:NumChans
        subplot(4,5,ch);
        plot(fF,PowSpect(:,rtms,ch),'Color',[1 .7 .7]); hold on;
        plot(fF,PowSpect(:,ltms,ch),'Color',[.7 .7 1]);
        plot(fF,mean(PowSpect(:,rtms,ch),2),'r'); hold on;
        plot(fF,mean(PowSpect(:,ltms,ch),2),'b');
    end
end
disp(['Time to plot (2) ' num2str(toc)]);

%% Classification
Ctypes = {'SW','L','SS','LO'};
% Findex = fF>=Frange(1) & fF<=Frange(2);
Frange = [6 30];
Findex = ismember(fF ,Frange(1):2:Frange(2));
Pfeat = [];
class =[];
PfeatAvg = [];
classAvg = [];
trlct = [];

if BuildClassifier == 1
    
    if Cfier == 3
        TwoClass = [2 1];
    else
        TwoClass = [1 -1];
    end
    
    
    pfolder = ls([bcifold '\parms\MuTask*']);
    pfolder = [bcifold '\parms\' pfolder '\mutask_' capt '.prm'];
    
    fpm = fopen(pfolder,'r');
    dmmy = 0;
    if fpm == -1
        disp('Parameter file does not exist!');
        fpm = fopen([bcifold  '\parms\' pfolder '\dummyfile.prm'],'r'); dmmy = 1;
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
        ImpInd(ll) = ~isempty(strfind(fpm2{1}{ll},'AcquisitionMode'));
        STVYInd(ll) = ~isempty(strfind(fpm2{1}{ll},'Sensitivity'));
        %         FIInd(ll) = ~isempty(strfind(fpm2{1}{ll},'FFTInputChannels'));
        NOInd(ll) = ~isempty(strfind(fpm2{1}{ll},'NormalizerOffsets'));
        NGInd(ll) = ~isempty(strfind(fpm2{1}{ll},'NormalizerGains'));
    end
    LinInd = find(LinInd==1);
    SNInd = find(SNInd==1); SNInd = SNInd(1);  %Needed to add for 4903BCI
    SEInd = find(SEInd==1); SEInd = SEInd(1);  %Needed to add for 4903BCI
    SFInd = find(SFInd==1);
    %     FIInd = find(FIInd==1);
    NOInd = find(NOInd==1);
    NGInd = find(NGInd==1);
    AccEvInd = find(AccEvInd==1);
    MinEvInd = find(MinEvInd==1);
    ImpInd = find(ImpInd==1);
    SBegInd = strfind(fpm2{1}{SNInd},'='); SEndInd = strfind(fpm2{1}{SNInd},'Name %');
    SEBegInd = strfind(fpm2{1}{SEInd},'='); SEEndInd = strfind(fpm2{1}{SEInd},'% %');
    SFBegInd = strfind(fpm2{1}{SFInd},'='); SFEndInd = strfind(fpm2{1}{SFInd},'//');
    BegInd = strfind(fpm2{1}{LinInd},'}'); EndInd = strfind(fpm2{1}{LinInd},'//');
    BegInd2 = strfind(fpm2{1}{LinInd},'= '); EndInd2 = strfind(fpm2{1}{LinInd},' {');
    %     FIBegInd = strfind(fpm2{1}{FIInd},'='); FIEndInd = strfind(fpm2{1}{FIInd},'//');
    NOBegInd = strfind(fpm2{1}{NOInd},'='); NOEndInd = strfind(fpm2{1}{NOInd},'% %');
    NGBegInd = strfind(fpm2{1}{NGInd},'='); NGEndInd = strfind(fpm2{1}{NGInd},'% %');
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
    
    
    TwoCues = {L_trialstart, L_trialend; R_trialstart R_trialend};
    
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
            ff = tT>(TwoCues{cc,1}(tt)/fs+.5) & tT<(TwoCues{cc,2}(tt)/fs-1.5);
            sum(ff)
            for ch = 1:NumChans
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
    
    if Figures_On==1
        figure('Name','Classifier features');
        if Avv == 0
            for i = 1:NumChans
                nrows = ceil(NumChans/4);
                ncols = ceil(NumChans/nrows);
                subplot(nrows,ncols,i); errorbar(fF(Findex),mean(Pfeat(class==TwoClass(2),...
                    (i-1)*sum(Findex)+1:i*sum(Findex))),...
                    std(Pfeat(class==TwoClass(2), (i-1)*sum(Findex)+1:i*sum(Findex))),'r');
                hold on
                errorbar(fF(Findex),mean(Pfeat(class==TwoClass(1),...
                    (i-1)*sum(Findex)+1:i*sum(Findex))),...
                    std(Pfeat(class==TwoClass(1), (i-1)*sum(Findex)+1:i*sum(Findex))),'b');
            end
        else
            for i = 1:NumChans
                nrows = ceil(NumChans/4);
                ncols = ceil(NumChans/nrows);
                subplot(nrows,ncols,i); errorbar(fF(Findex),mean(PfeatAvg(classAvg==TwoClass(2),...
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
    if Avv == 1
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
    
    if Figures_On==1
        if NumChans<=6
            figure('Name','Classifier features with outliers highlighted');
            if Avv==0
                for i = 1:NumChans
                    %Made a plotting change to visualize the windowed spectra that
                    %are removed.  Did not make this plotting change for the AvgF =
                    %1 case.
                    nrows = ceil(NumChans/4);
                    ncols = ceil(NumChans/nrows);
                    
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
                            (i-1)*sum(Findex)+1:i*sum(Findex)),'r'); end;
                    hold on;
                    if sum(ismember(temp1,tempK))~=0
                        plot(fF(Findex),Pfeat(temp1(ismember(temp1,tempK)),...
                            (i-1)*sum(Findex)+1:i*sum(Findex)),'b'); end;
                    
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
                for i = 1:NumChans
                    nrows = ceil(NumChans/4);
                    ncols = ceil(NumChans/nrows);
                    
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
    end
    
    
    %Remove high power (artifactual) features
    if Avv == 1
        PfeatAvg = PfeatAvg(~OB2,:);
        classAvg = classAvg(~OB2);
    else
        Pfeat = Pfeat(~OB2,:);
        class = class(~OB2);
    end
    trlct = trlct(~OB2);
    
    
    %DONT USE AVV ANYMORE
    for s2k = 1:5
        epochskip = 6-s2k;
        traind = Pfeat(1:epochskip:end,:);
        classlbl = class(1:epochskip:end);
        trlcc = trlct(1:epochskip:end);
        if Cfier == 1
            %Stepwise Linear Regression
            
            penter_v = .005; premove_v = .01;
            while penter_v<.9
                penter_v = penter_v*2;
                premove_v = premove_v*2;
                if premove_v>1
                    penter_v=.95;
                    premove_v = .99;
                end
                [B, SE, PVAL, INMODEL, STATS] = stepwisefit(traind,classlbl,...
                    'penter',.001,'premove',.01,'display','off');
                Classweight = -B(INMODEL==1);
                SigLocs = find(INMODEL==1);
            end
            
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
            
        elseif Cfier == 3
            
        elseif Cfier == 4 %LASSO
            
        elseif Cfier == 5 %Regularized LDA
            
            Mdl = fitcdiscr(traind,classlbl,...
                'SaveMemory','on');
            [err,gamma,delta,numpred] = cvshrink(Mdl,...
                'NumGamma',9,'NumDelta',9,'Verbose',1);
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
            
            
            
        end
        
        %Construct Classifier Matrix
        CM{s2k} = zeros(length(SigLocs),4);
        if isempty(CM{s2k})
            disp('No significant features!');
            CM{s2k} = [];
            skiptest = 1;
        else
            CM{s2k}(:,1) = (ceil(SigLocs./sum(Findex))'); %Input Channel
            Ctms = mod(SigLocs,sum(Findex));  Ctms(Ctms==0)=sum(Findex);
            Ffd = find(Findex==1);
            freqs = fF(Ffd);
            %         CM(:,2) = mat2cell(fF(Ffd(Ctms)),ones(size(SigLocs,2),1),1); %Input element (sample #)
            CM{s2k}(:,2) = fF(Ffd(Ctms));
            CM{s2k}(:,3) = 1; %Output channel
            CM{s2k}(:,4) = Classweight; %Weight
            % Changed from CM(:,4) = Classweight(SigLocs);  to reflect changes
            % to SW classifier that limits features to a specified number.
            
            
            
            
        end
        waitbar(.6+.3*s2k/5,wtbrr);
        
        
        %Looking at the expected classification for individual segments
        ffrt = fF(Findex);
        for ii = 1:size(traind,1)
            curs{s2k}(ii) = 0;
            for jj = 1:size(CM{s2k},1)
                chnn = CM{s2k}(jj,1);
                chind = (chnn-1)*length(ffrt)+1:chnn*length(ffrt);
                ffrx = ffrt==CM{s2k}(jj,2);
                curs{s2k}(ii) = curs{s2k}(ii) + CM{s2k}(jj,4)*traind(ii,chind(ffrx));
            end
        end
        for jj = unique(trlcc)'
            accum{s2k}(jj) = mean(curs{s2k}(trlcc==jj));
        end
        
    end
    temp = cellfun(@(x) ~isempty(x),CM);
    bigCM = max(find(temp));
    cm_full = CM{bigCM};
    

        %Looking at the expected classification for trial averages
 
    for ii = 1:size(PfeatAvg,1)
        cursA(ii) = 0;
        for jj = 1:size(CM{end},1)
            chnn = CM{end}(jj,1);
            ffr = (chnn-1)*sum(Findex)+1:chnn*sum(Findex);
            ffrx = ffrt==CM{end}(jj,2);
            cursA(ii) = cursA(ii) + CM{end}(jj,4)*PfeatAvg(ii,ffr(ffrx));
        end
    end
    


    %What is optimal decision boundary?
    %Set the gain and the offset of the Normalizer
    
    lcurs = curs{end}(class==TwoClass(1));
    rcurs = curs{end}(class==TwoClass(2));
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
    NormalizerGain = 1/std(curs{end});
    
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
        hist(curs{end}(class==TwoClass(1)),100);
        hold on;
        hist(curs{end}(class==TwoClass(2)),100);
        h = findobj(gca,'Type','patch');
        set(h(1),'FaceColor','r','EdgeColor','w','FaceAlpha',.5);
        set(h(2),'FaceColor','b','EdgeColor','w','FaceAlpha',.5);
        line([NormalizerOffset NormalizerOffset],[0 100],'Color','k','LineWidth',3)
    end
    
    
    
  
    
    %Instead of the line on the histogram plot, do an ROC curve
    %If we consider "truth" as being a left trial
    if ~isempty(Offrange)
    k=1;
    for thr = Offrange
        tp(k) = sum(rcurs>thr); fn(k) = sum(rcurs<thr);
        tn(k) = sum(lcurs<thr); fp(k) = sum(lcurs>thr); k = k+1;
    end
    tpr = tp./(tp+fn);
    fpr = fp./(fp+tn);
    auckinda = tpr.*(1-fpr);
    [~, emaxi] = max(auckinda);
    Offrange(emaxi);
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
    
    
    
    
    
    
        %%Append Hz and convert to savable form
    CM2 = mat2cell(CM{end},ones(size(CM{end},1),1),ones(size(CM{end},2),1));
    CM2(:,2) = cellfun(@num2str,CM2(:,2),'UniformOutput',false);
    for ii = 1:size(CM2,1);
        CM2{ii,2} = [CM2{ii,2} 'Hz'];
    end
    fid = fopen([bcifold '\data\' Name Sessions{max_sess_i} '\'...
        Name 'S' Sessions{max_sess_i} 'R' NewRunName '_' Ctypes{Cfier} 'Classifier.txt'],'wt');
    if fid==-1
        fid = fopen(['P:\ALS Proj Data\' Name Sessions{max_sess_i} '\'...
            Name 'S' Sessions{max_sess_i} 'R' NewRunName '_' Ctypes{Cfier} 'Classifier.txt']);
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
    fpm2{1}{LinInd} = [fpm2{1}{LinInd}(1:BegInd2+1) num2str(size(CM2,1)) ...
        fpm2{1}{LinInd}(EndInd2:BegInd+1) prmcls...
        fpm2{1}{LinInd}(EndInd-1:end)];
    fpm2{1}{LinInd} = horzcat(fpm2{1}{LinInd}{:});
%     %Change FFT Input Channels
%     fpm2{1}{FIInd} = [fpm2{1}{FIInd}(1:FIBegInd+1) num2str(size(SpatFilt,1)) ...
%         ' ' num2str(1:size(SpatFilt,1)) fpm2{1}{FIInd}(FIEndInd-1:end)];
    %Change Normalizer Offsets
    fpm2{1}{NOInd} = [fpm2{1}{NOInd}(1:NOBegInd+1) ...
        num2str([2 NormalizerOffset 0 0]) fpm2{1}{NOInd}(NOEndInd-1:end)];
    %Change Normalizer Gains
    fpm2{1}{NGInd} = [fpm2{1}{NGInd}(1:NGBegInd+1) ...
        num2str([2 NormalizerGain 1 0]) fpm2{1}{NGInd}(NGEndInd-1:end)];
    
    prmfile = [bcifold '\data\' Name Sessions{max_sess_i} '\mutask_'...
        Name 'S' Sessions{max_sess_i} 'R' NewRunName '.prm'];
    dlmcell(prmfile,fpm2{1});

    
else %Testrun - do not build classifier, prompt for a previously built one
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
    
    % Name = 'AndrewGrid-gUSB';
    % CSession = '002';
    % CRun = '02';
    % CM = dlmread(['C:\Documents and Settings\amg5106\My Documents\'...
    %         'Year2\Research\BCI2000\BCI2000sourcetree\data\AndrewGrid-gUSB002\'...
    %         'AndrewGrid-gUSBS002R01_SWClassifier.txt'],'\t');
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
end


%     CM(size(CM,1)+1,:) = [NormalizerGain 0 0 0];
%     CM(size(CM,1)+1,:) = [NormalizerOffset 0 0 0];

%% New Feedback with new classifier
%Cursor speed (DasherSpeller.cpp) is based on the Feedback interval length
if Figures_On == 1
    figure
    title('Left is 1, Right is 2');
    plot(1000*Feedback,'r');hold on;
    plot(999*TargCode,'b');
    plot(Cursor,'k');
%     if BuildClassifier == 0
        Classifiers = {ClassifierUsed, CM{end}};
        Normvalues = {[c.NormalizerGains.NumericValue(1) c.NormalizerOffsets.NumericValue(1)]...
            [NormalizerGain NormalizerOffset]};
%     elseif BuildClassifier == 1
%         Classifiers = {CM{end}};
%         Normvalues = {[NormalizerGain NormalizerOffset]}
%     end;
    
    colors = {'g','m'};
    for cv = 1:length(Classifiers)
        disp(['Classifier' num2str(cv)]);
        clsfr = Classifiers{cv};
        NormV = Normvalues{cv};
        CursorSpeed = 100 / (c.FeedbackDuration.NumericValue*(fs/sbs)) / 2;
        StartingPos = double([50]);
        TargetPos = c.Targets.Value(:,1:3);
        LTv = str2double(TargetPos{1,1}); RTv = str2double(TargetPos{2,1});
        CursorO = double(StartingPos(1)*ones(1,size(PowSpect,2)));
        tindex = wndw;
        targethit = 0;
        % for tt =wndw/8+1:95
        for tt = ceil(wndw/sbs+.5):size(PowSpect,2)
            %Is feedback running
            if Feedback(tindex) == 1
                %Did feedback just began, reset Cursor to Starting Position
                if Feedback(tindex-sbs) == 0
                    CursorO(:,tt) = StartingPos;
                    targethit = 0;
                else
                    %Increment Cursor
                    if targethit == 0
                        cursinc= double(zeros(3,1));
                        cursss(tt) = 0;
                        for i = 1:size(clsfr,1)
                            controlsignal = double(clsfr(i,4)*...
                                PowSpect(fF==clsfr(i,2),tt-wndw/sbs,clsfr(i,1)));
                            cursinc(clsfr(i,3)) = cursinc(clsfr(i,3))...
                                + controlsignal;
                            %                         cursss(tt) = cursss(tt) + clsfr(i,4)*...
                            %                             PowSpect(find(fF==clsfr(i,2)),tt-wndw/8,clsfr(i,1));
                        end
                        %Normalizer
                        cursinc = NormV(1)*(cursinc(1)-NormV(2));
                        
                        CursorO(:,tt) = CursorO(:,tt-1)+CursorSpeed*cursinc;
                        
                 
                        
                        %Check if a target has been hit
                            if TargCode(tindex) == 1 %right target is displayed
                                if CursorO(tt) >= RTv
                                    disp([num2str(tindex) 'in R_t'])
                                    targethit = 1;
                                end
                            elseif TargCode(tindex) == 2 %left/cross target
                                if CursorO(tt) <= LTv
                                    disp([num2str(tindex) 'in L_t'])
                                    targethit = 1;
                                end
                            else
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
            tindex = tindex+sbs;           
            
        end
        
        coordToState = 40.95;  %This is from the DoFeedback func. of DasherSpeller.cpp
        CursorO = fix(CursorO.*coordToState);
        %Only interested in the first index - the Y position of the cursor
        CursorO = repmat(CursorO,sbs,1);
        Output4 = reshape(CursorO,1,size(CursorO,1)*size(CursorO,2));
        
        
        plot(Output4,'Color',colors{cv})
    end
    
    
    
end



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

guess = zeros(length(unique(trlct)),s2k);
guess(cell2mat(cellfun(@(x) x>NormalizerOffset,accum,'UniformOutput',false)))=-1;
guess(cell2mat(cellfun(@(x) x<=NormalizerOffset,accum,'UniformOutput',false)))=1;

acc_fromCfier = sum(guess==repmat(classAvg,1,s2k))/length(classAvg);




sprintf(['Accuracies may not be the same because acc_fromCfier is \n'...
    'calculated from a fixed (3s) window, where the one calculated from\n'...
    'the plot is may be shorter if the target was reached sooner.'])


close(wtbrr)
end

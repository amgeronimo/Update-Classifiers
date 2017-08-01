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
% 12/15/15 - Modified v3 to allow for runs from multiple sessions to be
% used in classifier
%
% ~4/1/16 created BuildClassifier_teleBCI.m
%
% 4/20/16 - Changed code to allow for use of accumulate evidence.  Deals
% with the ability to aggregate flashes from multiple rounds of flashes
% 6/1/16 - teleBCI2 takes the location of the data as an input


function [CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD,ClassLabel, ChanLabel, SC, SR] = ...
    BuildClassifier_teleBCI2(Name,Sessions,Runs,BuildClassifier,Cfier,...
    Batch,RerefVal,Art,DoCSP,CSession,CRun, Figures_On, bcifold)

if Figures_On
    close all;
else
end
    wtbrr = waitbar(0,'Creating Classifier...');


disp('MAKE SURE ALIGN CHANNELS IS 0')
disp('When recording, data is saved after being filtered, ')
disp('no other processing has gone into saved data');
%Load Data and  Log File

cdp = cd;
[~, sl] = regexp(cdp,'\Users\');
se = regexp(cdp,'\');
tmp = find(sl==se);
pcusr = cdp(sl+1:se(tmp+1)-1);
if strcmp(pcusr,'ageronimo')
    %     bcifold = '\Documents\BCI2000_305_VS2012';
    %     bcifold = '\Documents\BCI2000_4903';
    bcifold = '\Documents\BCI2000_5300';
else
    bcifold = '\Documents\BCI2000_5341';
end


Avv = 0; %Classify from individual trials(0) or averages of stimulus codes(1)
% Figures_On = 0;


if Batch == 0
    RerefVal  = 0;  %0 = no reref, 1 = CAR, 2 = Lap, 3 = Bipolar (Neither Lap nor Bip work well)
    Art=0; %Artifact Data (0 - nothing, 1 - remove regions of artifact, 2 - regression)
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
     NameBase = ['C:\Documents and Settings\' pcusr bcifold '\data\' Name CSession...
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
NewRunName = [];
TTSpell = [];
GazeX = [];
GazeY = [];
SequencePhase = [];

kk = 1;
for ss = 1:length(Sessions)
    Session = Sessions{ss};
    for rr = 1:length(Runs{ss})
        %     [a b c d] = load_bcidat(['C:\Documents and Settings\amg5106\My Documents\'...
        %         'Year2\Research\BCI2000\BCI2000sourcetree\data\' ...
        %         Name Session '\' Name 'S' Session 'R' Run{rr} '.dat'],'-calibrated');
        try %try first to find the data on the computer
            [a b c d] = load_bcidat(['C:\Users\' pcusr bcifold '\data\' ...
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
            capt = '_emotiv';
            ChannelNames = cell(14,1);
            EEGloc = 1:14; EOGloc = [];
        elseif size(a,2)  == 22 %standard EEG cap (2 gusb)
            capt = '_22_3EOG_2amp';
            ChannelNames = {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8',...
                'P7','P3','Pz','P4','P8','O1','O2'};
            EEGloc = 1:19; EOGloc = 20:22;
        elseif size(a,2)  == 16 % 16 chan gNautilus
            capt = '_gNautilus16';
            ChannelNames = {'Fp1','Fp2','F3','Fz','F4','T7','C3','Cz','C4','T8',...
                'P3','Pz','P4','PO7','PO8','Oz'};
            EEGloc = 1:16; EOGloc = [];
        elseif size(a,2)  == 8 % 8 chan gNautilus
            capt = '_gNautilus8';
            ChannelNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};
            EEGloc = 1:8; EOGloc = [];
        elseif size(a,2)  == 3 
            capt = '_3ch';
            ChannelNames = cell(3,1);
        else
            disp('Undefined electrode configuration!')
            return
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
        
        
        
        fs = c.SamplingRate.NumericValue;
        %Trim Data and StimulusCodes
        trimMax = round((c.PreRunDuration.NumericValue-.5)*fs);
        Data = [Data;  a(trimMax:end,:)];
        StimulusCode = [StimulusCode; b.StimulusCode(trimMax:end)];
        StimulusType = [StimulusType; b.StimulusType(trimMax:end)];
        SourceTime = [SourceTime; b.SourceTime(trimMax:end)];
        StimulusTime = [StimulusTime; b.StimulusTime(trimMax:end)];
        TTSpell = [TTSpell cell2mat(c.TextToSpell.Value)];
        TTSpell_cell{kk} = cell2mat(c.TextToSpell.Value);
        if isfield(b,'EyetrackerLeftEyeGazeX')
            GazeX = [GazeX; b.EyetrackerLeftEyeGazeX(trimMax:end)];
            GazeY = [GazeY; b.EyetrackerLeftEyeGazeY(trimMax:end)];
        else
            GazeX = [GazeX; NaN(size(a(trimMax:end,:),1),1)];
            GazeY = [GazeY; NaN(size(a(trimMax:end,:),1),1)];
        end
        SequencePhase = [SequencePhase; b.PhaseInSequence(trimMax:end)];
        NewRunName = strcat(NewRunName,Runs{ss}{rr})
        if max(StimulusCode)==16  %Is this the multiflash grid (this) or single flash (next)?
            %                 fid = fopen(['C:\Documents and Settings\amg5106\My Documents\'...
            %                     'Year2\Research\BCI2000\BCI2000sourcetree\data\' Name Session '\'...
            %                     Name 'S' Session 'R' num2str(str2num(Run{rr})) '+_mylogfile.txt'],'r');
            
            
            fid = fopen(['C:\Documents and Settings\' pcusr bcifold '\data\' Name Session '\'...
                Name 'S' Session 'R' num2str(str2num(Runs{ss}{rr})) '+_mylogfile.txt'],'r');
            if fid==-1
                fid = fopen(['P:\ALS Proj Data\' Name Session '\'...
                    Name 'S' Session 'R' num2str(str2num(Runs{ss}{rr})) '+_mylogfile.txt'],'r');
            end
            
            itt=textscan(fid,'%s','delimiter','\n');
            InputText = [InputText; itt{1}];
            Speller = 1; %multiflash grid
            fclose(fid);
        elseif max(StimulusCode)==12 %The unmodified (RC) P300 speller
            Speller = 3;
        end
        
        %Also check for an RSVP speller, which can have any number of stimulus
        %codes, but will always have NumMatrixRow = 1. max(StimulusCode) will
        %equal 2*NumColumns
        if c.NumMatrixRows.NumericValue == 1 %This is the RSVP speller
            Speller = 2;
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
        NumTrials{kk} = fix((sum(b.StimulusCode>=1)/StateDuration)/(NumSequences{rr}*NumStimCodes));
        
        %Load Spatial Filter
        SpatFiltUsed{kk} = c.SpatialFilter.NumericValue;
        
        %Load Classifier
        ClassifierUsed{kk} = c.Classifier.NumericValue;
        ClassifierUsed{kk}(:,4) = cellfun(@str2num,c.Classifier.Value(:,4));
        
        kk = kk+1;
    end
end
    Data = Data';
    range = (round(.4*fs):round(fs)); %Use this data range for classifiation
    
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
    
    %Do the channel labels match the number of channels in Data?
    if length(EOGloc)+length(EEGloc)~=size(Data,1)
        disp('Channel labels incorrect'); return
    end
    
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
    
    %Are number of sequences the same
    if sum(cell2mat(NumSequences))/(kk-1) ~= NumSequences{1}
        disp('Not the same number of sequences in each run, code will not run correctly');
    end
    NumSequences = NumSequences{1};
    if SequencestoRemove >= NumSequences
        disp('Removed too many sequences');
    end
    
    %Number of rows, columns, stimuli
    NR = c.NumMatrixRows.NumericValue;
    NC = c.NumMatrixColumns.NumericValue;
    NS = size(c.TargetDefinitions.NumericValue,1);
    
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
        plot(StimulusCode(1:c.SampleBlockSize.NumericValue:end),'g');
        plot(-BadTime*10,'k')
    end
    %figure
    %plot(Data(:,1)*1000); hold on; plot(b.SourceTime,'r');
    
    
    %% Find superlarge (>300 mV) artifacts
    
    %If not using the contents of this cell, set GoodData to ones;
    GoodData = ones(size(Data,2),1);
    
    
    BadData = [];
    for i = 1:size(Data,1)
        %Find large artifacts
        BadData = [BadData find(Data(i,:)>300|Data(i,:)<-300)];
    end
    %Define a range 1 second before and 3 seconds after all of the found
    %artifacts
    %Commented this out int the haste of recording.  Was giving an error for P31_CrossCBS002R02
    % BadData = unique(BadData)
    % BadData = repmat(BadData,4*fs,1)+repmat((-1*fs+1:3*fs)',1,length(BadData));
    % BadData = unique(BadData);
    GoodData = ~ismember(1:length(Data),BadData);
    
    %% Plot raw data and trim bad data
    
    if Figures_On == 1
        figure('Name','Data with artifacts marked')
        for i = 1:size(Data,1)
            if size(Data,1) == 22
                eleclocs = [7 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 29 2 3 4];
                tr(i) = subplot(6,5,eleclocs(i)); plot(Data(i,:));
                hold on;
                plot(100*(1-GoodData),'r');
            else
                tr(i) = subplot(6,5,i); plot(Data(i,:));
                hold on;
                plot(100*(1-GoodData),'r');
            end
        end
        linkaxes(tr);
    end
    
    GoodData = logical(GoodData);
    
    
    %Not trimming bad data -- too complicated to change all the remaining code
    %to work with non-regular stimulus codes (much is based on having exactly
    %16 codes of 10 repetitions each.  ***Performing Epoch rejection later***
    %
    % Data = Data(:,GoodData);
    % StimulusCode = StimulusCode(GoodData);
    % StimulusType =  StimulusType(GoodData);
    % SourceTime = SourceTime(GoodData);
    % StimulusTime = StimulusTime(GoodData);
    %
    % figure
    % plot(Data(1,:),'r'); hold on; plot(StimulusCode,'g','LineWidth',2);
    % GDi = find(GoodData-circshift(GoodData,[0 1])==-1)
    % line([GDi GDi],[-100 100])
    
    %% Artifact Correction
    
    if Art == 1 %Remove data displaying artifacts
        try
            ArtInt = find(Data(EOGloc(1),:)>70&circshift(Data(EOGloc(1),:),[0 1])<70);
        catch
            disp('Cannot remove data because no EOG channels were recorded');
            return
        end
        Arange = fix(-.1*fs:1*fs);
        ArtInt = repmat(ArtInt,length(Arange),1)+repmat(Arange,length(ArtInt),1)';
        %     ArtInt = ArtInt';
        ArtInt = ArtInt(:);
        ArtInt = ArtInt(ArtInt>0 & ArtInt<size(Data,2));
        ArtLoc = zeros(size(Data,2),1);
        ArtLoc(ArtInt)=1;
        
        if Figures_On == 1
            figure('Name','EOG channel 1, locations of artifact, and stimulus code')
            plot(Data(EOGloc(1),:)); hold on; plot(50*ArtLoc,'r'); plot(StimulusCode,'g');
            %     figure
            %     tt(1) = subplot(211); plot(Data(EEGloc,:)');
            %     tt(2) = subplot(212); plot(Data(EOGloc,:)');
            %     linkaxes(tt)
        end
        
        %This is the data with the Artifacts removed
        Data = Data(:,~ArtLoc);
        %Do the same for other states as well (Taking care of triggers below)
        StimulusType = StimulusType(~ArtLoc);
        StimulusTime = StimulusTime(~ArtLoc);
        SourceTime = SourceTime(~ArtLoc);
        StimulusCode = StimulusCode(~ArtLoc);
        
        if Figures_On == 1
            figure('Name','Channel 14 and stimulus code')
            plot(Data(14,:)); hold on; plot(StimulusCode,'g');
        end
        %Now shift the triggers back
        
        
        ArtWeights = eye(size(Data,1));
        ArtWeights2 = eye(length(EEGloc),size(Data,1));
        NumChans = length(EEGloc);
        
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
        
    else
        ArtWeights = eye(size(Data,1));
        ArtWeights2 = eye(length(EEGloc),size(Data,1));
        NumChans = length(EEGloc);
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
    
    NumChans;
    
    disp(['Done with Art ' num2str(toc)]);
    
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
            ChosenLind2 = find(ChosenLind)+8; %this is the distance between the selected letter and the repeat of the next goal letter.
            trialkeep = ~ismember(GoalLind2,ChosenLind2); 
            trialkeep(end)=0;
            GoalLindF = [GoalLindF (GoalLind2(trialkeep))];
            GoalLetter = InputText(GoalLind2(trialkeep));
            GoalLetter = cell2mat(GoalLetter);
            GoalTrial = str2num(GoalLetter(:,18:19));
            GoalTrialF = [GoalTrialF; (GoalTrial+max(GoalTrialF))];%This 
            %vector is used to reconcile the number of trials and the number of letters spelled
            GoalLetter = TTSpell_cell{rr}(GoalTrial)';
            GoalLetterF = [GoalLetterF; GoalLetter];
            ChosenLetter = cell2mat(InputText(ChosenLind(1:end-1)));
            ChosenLetter = ChosenLetter(:,18);
            ChosenLetterF = [ChosenLetterF; ChosenLetter];
            Begind = find(Begin==1);
            SCind = find(SCind==1)'; FSCi = find(SCind < Begind);
            SCindF = [SCindF; SCind(FSCi(end):end-1)];
            SRind = find(SRind==1)'; FSRi = find(SRind < Begind);
            SRindF = [SRindF; SRind(FSRi(end):end-1)];
            INITmind = find(INITmind==1)';
            INITmindF = [INITmindF; INITmind(1:end-1)];
            NEXTmind = find(NEXTmind==1)';  omitN = NumSequences:...
                NumSequences:NumTrials{rr}*NumSequences;
            NEXTmind(omitN) = [];
            NEXTmindF = [NEXTmindF; NEXTmind];
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
    
    
        
        NumTrials = sum(cell2mat(NumTrials));
        NEXTm = cell(NumSequences,length(GoalLetterF));
        for i = 1:length(GoalLetterF)
            tempSC = InputText{SCindF(GoalTrialF(i))}(4:end);
            tempSC = regexp(tempSC,' ','split')'; tempSC(cellfun(@isempty,tempSC)) = [];
            SC{i} = tempSC;
            tempSR = InputText{SRindF(GoalTrialF(i))}(4:end);
            tempSR = regexp(tempSR,' ','split')'; tempSR(cellfun(@isempty,tempSR)) = [];
            SR{i} = tempSR;
            tempINITm = InputText{INITmindF(i)}(16:end);
            tempINITm = regexp(tempINITm,' ','split')'; tempINITm(cellfun(@isempty,tempINITm)) = [];
            INITm{i} = tempINITm;
            for j = 1:NumSequences-1
                tempNEXTm = InputText{NEXTmindF((NumSequences-1)*(i-1)+j)}(16:end);
                tempNEXTm = regexp(tempNEXTm,' ','split')'; tempNEXTm(cellfun(@isempty,tempNEXTm)) = [];
                NEXTm{j,i} = tempNEXTm;
            end
        end
        if isempty(NEXTm)
            mSeq = INITm;
        else
            mSeq = cat(1,INITm,NEXTm);
        end
        %%%%%%%%%%%
        
    elseif Speller == 2
        NumTrials = sum(cell2mat(NumTrials));
        for i = 1:NumTrials
            SC{i} = 1:c.NumMatrixColumns.NumericValue;
            SR{i} = c.NumMatrixColumns.NumericValue+1:2*c.NumMatrixColumns.NumericValue;
        end
        
    elseif Speller == 3
        NumTrials = sum(cell2mat(NumTrials));
        for i = 1:NumTrials
            SC{i} = repmat(7:12,1,NR);
            SR{i} = repmat(1:6,NC,1);
            SR{i} = SR{i}(:)';
        end
    end
    
    
    %Check that the markers saved in the log file are the same as those in the
    %data file
    tStimC = double(nonzeros(StimulusCode)); tStimC = tStimC(1:StateDuration:end);
    for i = 1:NumTrials
        for j = 1:NumSequences
            StimCode{j,i} = num2cell(tStimC(NumStimCodes*NumSequences*(i-1)+NumStimCodes*(j-1)+1:...
                NumStimCodes*NumSequences*(i-1)+NumStimCodes*(j)));
        end
    end
    
    if Speller == 1
        %StimCode should contain the same values as mSeq
        ies=0;
        k = 1;
        for i = 1:NumTrials
            for j = 1:NumSequences
                ies = isequal(cell2mat(StimCode{j,i}),str2double(mSeq{j,i}))+ies;
                isgood(k) = isequal(cell2mat(StimCode{j,i}),str2double(mSeq{j,i}));
                k=k+1;
            end
        end
        if ies ~= NumTrials*NumSequences
            disp('SOMETHING IS WRONG')
        end
    end
    
    
    %If removing data(probably wont do) need to also remove the same portion of
    %StimCode

    NumSequences = NumSequences-SequencestoRemove;
    StimCode = StimCode(1:NumSequences,:);
    
    disp(['Done with myfile ' num2str(toc)]);
    waitbar(.1,wtbrr);
    
    
    %% Do CSP
    if DoCSP == 1
        if BuildClassifier == 1
            Data_C = (RefFilt*ArtWeights2*Data)';
            [W] = CSP_P300(Data_C,StimulusCode,StimulusType,...
                NumChans, NumTrials, NumStimCodes,range,StateDuration);
            [W2] = CSP_P300_ver2(Data_C,StimulusCode,StimulusType,...
                NumChans, NumTrials, NumStimCodes,range,StateDuration);
            
        elseif BuildClassifier == 0
        end
    elseif DoCSP == 0
        if BuildClassifier == 1
            W = eye(size(RefFilt,1));
            W2 = eye(size(RefFilt,1));
        else
        end
    end
    
    disp(['Done with CSP ' num2str(toc)]);
    
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
                Clocs = find(StimulusCode>0 & StimulusType == 1);
                Clocs = Clocs(1:StateDuration:end);
                Nlocs = find(StimulusCode>0 & StimulusType == 0);
                Nlocs = Nlocs(1:StateDuration:end);
                CCdata = zeros(length(range),size(Data,1),length(Clocs));
                NNdata = zeros(length(range),size(Data,1),length(Nlocs));
                CCdataA = zeros(length(range),size(DataA,1),length(Clocs));
                NNdataA = zeros(length(range),size(DataA,1),length(Nlocs));
                CCdataAF = zeros(length(range),size(DataAF,1),length(Clocs));
                NNdataAF = zeros(length(range),size(DataAF,1),length(Nlocs));
                CCdataAFC = zeros(length(range),size(DataAFC,1),length(Clocs));
                NNdataAFC = zeros(length(range),size(DataAFC,1),length(Nlocs));
                CCdataAFCT = zeros(length(range),size(DataAFCT,1),length(Clocs));
                NNdataAFCT = zeros(length(range),size(DataAFCT,1),length(Nlocs));
                toc
                for cl = 1:length(Clocs)
                    CCdata(:,:,cl) = Data(:,Clocs(cl)+range-1)';
                    CCdataA(:,:,cl) = DataA(:,Clocs(cl)+range-1)';
                    CCdataAF(:,:,cl) = DataAF(:,Clocs(cl)+range-1)';
                    CCdataAFC(:,:,cl) = DataAFC(:,Clocs(cl)+range-1)';
                    CCdataAFCT(:,:,cl) = DataAFCT(:,Clocs(cl)+range-1)';
                end
                toc
                for nl = 1:length(Nlocs)
                    NNdata(:,:,nl) = Data(:,Nlocs(nl)+range-1)';
                    NNdataA(:,:,nl) = DataA(:,Nlocs(nl)+range-1)';
                    NNdataAF(:,:,nl) = DataAF(:,Nlocs(nl)+range-1)';
                    NNdataAFC(:,:,nl) = DataAFC(:,Nlocs(nl)+range-1)';
                    NNdataAFCT(:,:,nl) = DataAFCT(:,Nlocs(nl)+range-1)';
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
                for np = 1:size(p3p,1)
                    cdta = p3p{np,1}; ndta = p3p{np,2};
                    if Figures_On == 1
                    figure('Name',p3pn{np})
                    for ch = 1:NumChans
                        subplot(4,5,ch);
%                         plot(range/fs,squeeze(cdta(:,ch,:)),'Color',[.7 .7 1]);
%                         hold on
%                         plot(range/fs,squeeze(ndta(:,ch,:)),'Color',[1 .7 .7]);
                        plot(range/fs,nanmean(cdta(:,ch,:),3),'Color','b');
                        hold on
                        plot(range/fs,nanmean(ndta(:,ch,:),3),'Color','r');
                    end
                    end
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
            dlmwrite(['C:\Documents and Settings\' pcusr bcifold '\data\' Name Session '\'...
                Name 'S' Session 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
        catch
            dlmwrite(['P:\ALS Proj Data\' Name Session '\'...
                Name 'S' Session 'R' NewRunName '_SpatialFilter.txt'],SpatFilt,'delimiter','\t')
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
    if ~isempty(GazeX)
        [gazeacc, gazevar, gazeinvalid] = P300_EyeTracking_inP300Classifier_v2(...
            GazeX,GazeY, Data(EOGloc,:), SC,SR, StimulusCode, ...
            SequencePhase, StimulusType, StateDuration, ...
            NumSequences+SequencestoRemove, Speller,...
            fs, c, Figures_On, 0);
    end
    disp(['Finished with eye tracking ' num2str(toc)]);
    waitbar(.25,wtbrr);
    %% Process and Organize Data
    
    %%%%
    %Associate Target Codes with Targets
    TargSymb = c.TargetDefinitions.Value(:,2)';
    SymbSymb = {'Yes','Rest','Wash','Toil','Call','NO','Doc','Nurse','Hus','Fam',...
        'Vis','Wheel','Walk','Pillow','Blanket','Clothes','Med','Pray','Clerg',...
        'Drink','Eat','Slip','Pap','Tv','Music','Book','Glass','Breat','Cold',...
        'Hot','Time','Date'};
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Organize Data into choice and not-choice targets
    
    
    if FreeSpell==0
        targetstim = double(StimulusCode(logical(StimulusType)));
        targetstim = targetstim(1:StateDuration:length(targetstim));
        targetstim = [targetstim(1:(NumSequences+SequencestoRemove)*2:length(targetstim))...
            targetstim(2:(NumSequences+SequencestoRemove)*2:length(targetstim))];
    else
        targetstim = repmat([1 2],NumTrials,1); %Put fake target stims in... just
        %so the code below runs and we group the data for classification.
        %We will ignore this when classifying later
    end
    
    
    
    %Initialize counters for each stimulus code
    Tind = ones(1,NumStimCodes);
    AllData = zeros(NumStimCodes*NumTrials*NumSequences*NumChans,length(range));
    AllData_c = zeros(NumStimCodes*NumTrials*NumSequences*NumChans,1);
    AllData_t = zeros(NumStimCodes*NumTrials*NumSequences*NumChans,1);
    ad_i = 1;
    for i = 1:NumTrials
        for j = 1:NumStimCodes %Loop through stimulus Codes
            TargLoc = find(StimulusCode == j); TargLoc = TargLoc(1:StateDuration:end);
            
            if ismember(j,targetstim(i,:))
                tmrk = 1;
            else
                tmrk = 0;
            end
            %Now take the locations of the first NumSeq of each of these vectors, and
            %increment a counter by NumSeq
            TrialTargLoc{i,j} = TargLoc(Tind(j):Tind(j)+NumSequences-1);
            Tind(j) = Tind(j)+NumSequences+SequencestoRemove;
            
            for ch = 1:NumChans %Loop through channels
                TrialTargData{i,j}{ch} = [];
                TrialTargData_i{i,j}{ch} = [];
                %             figure
                for s = 1:NumSequences %Loop through sequences per trial, extract data 0-.6 seconds after stimulus
                    
                    TrialTargData{i,j}{ch} = [TrialTargData{i,j}{ch}; DD(ch,TrialTargLoc{i,j}(s)+range-1)];
                    TrialTargData_i{i,j}{ch} = [TrialTargData_i{i,j}{ch}; ad_i];
                    AllData(ad_i,:) = DD(ch,TrialTargLoc{i,j}(s)+range-1);
                    AllData_c(ad_i) = ch;
                    AllData_t(ad_i) = tmrk;
                    ad_i = ad_i+1;
                    %            plot(StimulusCode); hold on;
                    %            plot(TrialTargLoc{i,j}(s)+range(1)-1:TrialTargLoc{i,j}(s)+range(end)-1,...
                    %                 DD(TrialTargLoc{i,j}(s)+range(1)-1:TrialTargLoc{i,j}(s)+range(end)-1,ch)'); hold on;
                end
            end
        end
    end
    
    disp(['Time for organizing data ' num2str(toc)]);
    waitbar(.3,wtbrr);
    % %% Old method (slower)
    % tic
    % oEEGData = [];
    % oChanLabel = [];
    % oClassLabel = [];
    % oAvgEEGData = [];
    % oAvgChanLabel = [];
    % oAvgClassLabel = [];
    % for ch = 1:NumChans %Loop through channels
    %     odataloc{ch} = [];
    %     for i = 1:size(StimCode,2) %Loop through Trials
    %
    %         for j = targetstim(i,:)
    %             oEEGData = [oEEGData; TrialTargData{i,j}{ch}];
    %             oAvgEEGData = [oAvgEEGData; mean(TrialTargData{i,j}{ch},1)];
    %             oAvgChanLabel = [oAvgChanLabel; ch];
    %             oAvgClassLabel = [oAvgClassLabel; 'C'];
    %             oChanLabel = [oChanLabel; repmat(ch,NumSequences,1)];
    %             oClassLabel = [oClassLabel; repmat('C',NumSequences,1)];
    %             odataloc{ch} = [odataloc{ch}; TrialTargLoc{i,j}];
    %         end
    %         for k = find(~ismember(1:NumStimCodes,targetstim(i,:))==1)
    %             oEEGData = [oEEGData; TrialTargData{i,k}{ch}];
    %             oAvgEEGData = [oAvgEEGData; mean(TrialTargData{i,k}{ch},1)];
    %             oAvgChanLabel = [oAvgChanLabel; ch];
    %             oAvgClassLabel = [oAvgClassLabel; 'N'];
    %             oChanLabel = [oChanLabel; repmat(ch,NumSequences,1)];
    %             oClassLabel = [oClassLabel; repmat('N',NumSequences,1)];
    %             odataloc{ch} = [odataloc{ch}; TrialTargLoc{i,k}];
    %         end
    %     end
    % end
    % disp(['Time to rear data (1)' num2str(toc)]);
    
    %Changed from 'N' and 'C' based to 2 and 1 based.
    
    % EEGData = zeros(NumSequences*NumTrials*NumStimCodes*NumChans,length(range));
    % ChanLabel = zeros(size(EEGData,1),1);
    % ClassLabel = repmat(char(0),size(EEGData,1),1);
    % AvgEEGData = zeros(NumTrials*NumStimCodes*max(ChanLabel),length(range));
    % AvgChanLabel = zeros(size(AvgEEGData,1),1);
    % AvgClassLabel = repmat(char(0),size(AvgEEGData,1),1);
    % for ch = 1:NumChans %Loop through channels
    %     dataloc{ch} = [];
    %     for i = 1:NumTrials %Loop through Trials
    %
    %         for j = 1:NumStimCodes
    %             index = (ch-1)*(NumTrials*NumStimCodes*NumSequences)+(i-1)*...
    %                 (NumStimCodes*NumSequences)+(j-1)*NumSequences+1:(ch-1)*...
    %                 (NumTrials*NumStimCodes*NumSequences)+(i-1)*(NumStimCodes*...
    %                 NumSequences)+(j)*NumSequences;
    %             avg_index = (ch-1)*(NumTrials*NumStimCodes)+(i-1)*...
    %                 (NumStimCodes)+(j-1)+1:(ch-1)*(NumTrials*NumStimCodes)+...
    %                 (i-1)*(NumStimCodes)+(j);
    %             EEGData(index,:) = TrialTargData{i,j}{ch};
    %             AvgEEGData(avg_index,:) = mean(TrialTargData{i,j}{ch},1);
    %
    %
    %             if ismember(j,targetstim(i,:)) %Choice targets = 1
    %              ChanLabel(index) = repmat(ch,NumSequences,1);
    %             ClassLabel(index) = repmat('C',NumSequences,1);
    %             AvgChanLabel(avg_index) = ch;
    %             AvgClassLabel(avg_index) = 'C';
    %             else
    %                ChanLabel(index) = repmat(ch,NumSequences,1);
    %             ClassLabel(index) = repmat('N',NumSequences,1);
    %             AvgChanLabel(avg_index) = ch;
    %             AvgClassLabel(avg_index) = 'N';
    %             end
    %
    %
    %             dataloc{ch}(index) = TrialTargLoc{i,j};
    %         end
    %
    %     end
    % end
    
    
    
    %Find indices of bad (unaveraged) epochs here, then average sequences, for
    %each channel separately and for choice and non-choice separately.
    BadEpoch = zeros(NumStimCodes*NumTrials*NumSequences*NumChans,1);
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
    
    EEGData = zeros(NumSequences*NumTrials*NumStimCodes*NumChans,length(range));
    EEGDataRR = zeros(NumSequences*NumTrials*NumStimCodes*NumChans,length(range));
    ChanLabel = zeros(size(EEGData,1),1);
    ClassLabel = zeros(size(EEGData,1),1);
    AvgChanLabel = zeros(NumTrials*NumStimCodes*NumChans,1);
    AvgClassLabel = zeros(NumTrials*NumStimCodes*NumChans,1);
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
                            NumSequences)+(j)*NumSequences;
                        EEGData(idx,:) = [TrialTargData{i,j}{ch}(good_trl,:); NaN(NumSequences-length(good_trl),length(range))];
                        EEGDataRR(idx,:) = TrialTargData{i,j}{ch};
                    
                    end
                    
                    if ns>length(good_trl)
                        epc_to_avg = good_trl;
                    else
                        epc_to_avg = good_trl(1:ns);
                    end
                    
                    avg_index = (ch-1)*(NumTrials*NumStimCodes)+(i-1)*...
                        (NumStimCodes)+(j-1)+1:(ch-1)*(NumTrials*NumStimCodes)+...
                        (i-1)*(NumStimCodes)+(j);
                    
                    
                    AvgEEGData{ns}(avg_index,:) = mean(TrialTargData{i,j}{ch}(epc_to_avg,:),1);
                    AvgEEGDataRR{ns}(avg_index,:) = mean(TrialTargData{i,j}{ch});
                    %             if ch==1 && i==1 && j==1
                    %             disp('Detrending P300 epochs!');
                    %             end
                    %             EEGData(index,:) = detrend(TrialTargData{i,j}{ch},0);
                    %             AvgEEGData(avg_index,:) = detrend(mean(TrialTargData{i,j}{ch},1),0);
                    %
                    
                    %             SScodes(index)=j;
                    if ns == NumSequences
                        if ismember(j,targetstim(i,:)) %Choice targets = 1
                            ChanLabel(idx) = repmat(ch,NumSequences,1);
                            ClassLabel(idx) = ones(NumSequences,1);
                            AvgChanLabel(avg_index) = ch;
                            AvgClassLabel(avg_index) = 1;
                        else
                            ChanLabel(idx) = repmat(ch,NumSequences,1);
                            ClassLabel(idx) = 2*ones(NumSequences,1);
                            AvgChanLabel(avg_index) = ch;
                            AvgClassLabel(avg_index) = 2;
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
    if Figures_On == 1
        figure
        imagesc(EEGData); colorbar
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
                    plot(range/fs,mean(AvgEEGData{ns}(AvgClassLabel==2&AvgChanLabel==ch,:)),'Color','r');
                    plot(range/fs,mean(AvgEEGData{ns}(AvgClassLabel==1&AvgChanLabel==ch,:)),'Color','b');
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
                for ch = 1:NumChans
                    subplot(4,5,ch);
                    plot(range/fs,mean(EEGDataRR(ClassLabel==2&ChanLabel==ch,:)),'Color',[1 .7 .7]);
                    hold on;
                    plot(range/fs,mean(EEGDataRR(ClassLabel==1&ChanLabel==ch,:)),'Color',[.7 .7 1]);
                    plot(range/fs,nanmean(EEGData(ClassLabel==2&ChanLabel==ch,:)),'Color','r');
                    plot(range/fs,nanmean(EEGData(ClassLabel==1&ChanLabel==ch,:)),'Color','b');
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
        
        TrainD = [];
        for ch = 1:max(ChanLabel)
            TrainD = [TrainD double(x_down(ChanLabel==ch,:))];
            inpch = [inpch ch*ones(1,size(x_down,2))];
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
                    (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes)),'Color',[1 .7 .7]); hold on;
                plot(LDAtimes,AvgTrainD{ns}(AvgClassLabel(AvgChanLabel==1)==1,...
                    (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes)),'Color',[.7 .7 1]);
                plot(LDAtimes,mean(AvgTrainD{ns}(AvgClassLabel(AvgChanLabel==1)==2,...
                    (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes))),'r'); hold on;
                plot(LDAtimes,mean(AvgTrainD{ns}(AvgClassLabel(AvgChanLabel==1)==1,...
                    (ch-1)*length(LDAtimes)+1:ch*length(LDAtimes))),'b');
            end
        elseif FreeSpell == 1
        end
    end
    disp(['Time to plot (2) ' num2str(toc)]);
    
    %% Classification
    ClassTypes = {'SW','LDA','SS','LO','rLDA'};
    if BuildClassifier == 1
        
        
        %Open the parameter file for writing
        if Speller == 3 % Grid
            fpm = fopen(['C:\Documents and Settings\' pcusr bcifold '\parms\CBSpeller\P3Speller_calib' capt '.prm'],'r');
        elseif Speller == 1 %Checkerboard
            if c.IconHighlightMode.NumericValue == 5
                %Check if using shuffled faces
                if isempty(strfind(c.OverlayImagePaths.Value{1},'_s'))
                    fpm = fopen(['C:\Documents and Settings\' pcusr bcifold '\parms\CBSpeller\FACE_CrossSpeller_calib' capt '.prm'],'r');
                else
                    fpm = fopen(['C:\Documents and Settings\' pcusr bcifold '\parms\CBSpeller\FACEshuffle_CrossSpeller_calib' capt '.prm'],'r');
                end
            else
                fpm = fopen(['C:\Documents and Settings\' pcusr bcifold '\parms\CBSpeller\FLSH_CrossSpeller_calib' capt '.prm'],'r');
            end
        elseif Speller == 2 %Single Flash (grid or RSVP)
            if c.NumMatrixColumns.NumericValue == 8
                fpm = fopen(['C:\Documents and Settings\' pcusr bcifold '\parms\RSVPSpeller\RSVP_P3_8targ_calib' capt '.prm'],'r');
            else
                fpm = fopen(['C:\Documents and Settings\' pcusr bcifold '\parms\RSVPSpeller\RSVP_P3_calib' capt '.prm'],'r');
            end
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
            TTS(ll) = ~isempty(strfind(fpm2{1}{ll},'TextToSpell'));
            ImpInd(ll) = ~isempty(strfind(fpm2{1}{ll},'AcquisitionMode'));
            IntInd(ll) =  ~isempty(strfind(fpm2{1}{ll},'InterpretMode'));
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
            S2K = repmat(S2K,NumStimCodes*NumTrials,1);
            traind = TrainD(S2K,:);
            rem_ind = sum(isnan(traind),2)>0;
            traind = traind(~rem_ind,:);
            size(traind)
            classlbl = ClassLabel(ChanLabel==1);
            classlbl = classlbl(S2K);
            classlbl = classlbl(~rem_ind);
            if Cfier == 1
                %Stepwise Linear Regression
                INMODEL = zeros(1,size(AvgTrainD{s2k},2));
                penter_v = .005; premove_v = .01;
                while sum(INMODEL)<4
                    penter_v = penter_v*2
                    premove_v = premove_v*2
                    if penter_v>1
                        penter_v=.95
                        premove_v = .99
                    end
                    if Avv == 0
                        [B, SE, PVAL, INMODEL] = stepwisefit(traind,...
                            classlbl,'penter',penter_v,'premove',premove_v,'display','off');
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
                    [class,err,post,logp,coeff] = classify(traind, traind, classlbl);
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
        bigCM = max(find(temp))
        cm_full = CM{bigCM};
        
        
        
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
        [ma mi] = max(psums)
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
                    for kk = jj+1:size(TrialChoice,2) %Loop through code 2
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
                        elseif Speller == 2 || Speller == 3
                            if ismember(trlchc(jj,1),SR{i})
                                mtch = find(SR{i}==trlchc(jj,1)&...
                                    SC{i}==trlchc(jj,2));
                            else
                                mtch = find(SC{i}==trlchc(jj,1)&...
                                    SR{i}==trlchc(jj,2));
                            end
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
            
            WORD_full(s2k,:) = TargSymb(TargChoice);
            if FreeSpell == 0
                try
%                 WORDD_full(s2k) = sum(cell2mat(WORD_full(s2k,:)) == TTSpell)/length(TTSpell)*100;
                WORDD_full(s2k) = sum(cell2mat(WORD_full(s2k,:)) == GoalLetterF')/length(GoalLetterF)*100

                catch
                    disp('not calculating accuracy correctly')
                    WORDD_full(s2k) = NaN;
                end
            end
        end
    end
    
    
    waitbar(.9,wtbrr);
    
    SequencesUsed = (1:NumSequences)';
    Guess = cell2mat(WORD_full);
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
    WORD_Acc = Accuracy(end);
    %% Save parameter files
    if BuildClassifier == 1
        try
            dlmwrite(['C:\Documents and Settings\' pcusr bcifold '\data\' Name Session '\'...
                Name 'S' Session 'R' NewRunName '_' ClassTypes{Cfier} 'Classifier.txt'],cm_full,'delimiter','\t')
        catch
            dlmwrite(['P:\ALS Proj Data\' Name Session '\'...
                Name 'S' Session 'R' NewRunName '_' ClassTypes{Cfier} 'Classifier.txt'],cm_full,'delimiter','\t')
        end
        
        skiptest = 0;
        
        
        %save in parameter file
        prmcls = num2str(reshape(cm_full',1,size(cm_full,1)*size(cm_full,2)));
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
        fpm2{1}{LinInd} = [fpm2{1}{LinInd}(1:BegInd2+1) num2str(size(cm_full,1)) ...
            fpm2{1}{LinInd}(EndInd2:BegInd+1) prmcls...
            fpm2{1}{LinInd}(EndInd-1:end)];
        %WordBank
        WordBank = {'JAWS','LINK','PARK','TEN9','C3PO','R2D2','CART','MIND',...
            'ALOE','PINK','BLUE','GREY','TAPE','CAPE','BURN','WIND','WILD',...
            'JUMP','PUMP','CAMP','MAPS','BARN','COAL','BOOT','TOLL','HUNT',...
            'FREE','COOL','MARK'};
        temp = randperm(length(WordBank));
        %Change the Word to be spelled
        fpm2{1}{TTSInd} = [fpm2{1}{TTSInd}(1:TTSBegInd+1) WordBank{temp(1)} ...
            fpm2{1}{TTSInd}(TTSEndInd-1:end)];
        %Change the letters to be displayed
        fpm2{1}{DispResInd} = [fpm2{1}{DispResInd}(1:DispResBegInd+1) '1 1 0 1 ' ...
            fpm2{1}{DispResInd}(DispResEndInd:end)];
        %Add Minimum Evidence
        if ~isempty(AccEvInd)
            if MinEvidence(bigCM) < 0
                disp('Target evidence not strong enough to set MinEvidence parameter')
                MinEvidence(bigCM) = 0;
            end
            prmme = [num2str(MinEvidence(bigCM)) ' 0 0 % '];
            fpm2{1}{MinEvInd} = [fpm2{1}{MinEvInd}(1:MBegInd+1) prmme...
                fpm2{1}{MinEvInd}(MEndInd-1:end)];
        end
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
        try
            dlmcell(['C:\Documents and Settings\' pcusr bcifold '\data\' Name Session '\CopyParams'...
                Name 'S' Session 'R' NewRunName '_' ClassTypes{Cfier} '.prm'],fpm2{1});
            dlmcell(['C:\Documents and Settings\' pcusr bcifold '\data\' Name Session '\FreeParams'...
                Name 'S' Session 'R' NewRunName '_' ClassTypes{Cfier} '.prm'],fpm2f{1});
        catch
            dlmcell(['P:\ALS Proj Data\' Name Session '\ParamFile'...
                Name 'S' Session 'R' NewRunName '_' ClassTypes{Cfier} '.prm'],fpm2{1});
            dlmcell(['P:\ALS Proj Data\' Name Session '\FreeParams'...
                Name 'S' Session 'R' NewRunName '_' ClassTypes{Cfier} '.prm'],fpm2f{1});
        end
    end
    
close(wtbrr)   
end

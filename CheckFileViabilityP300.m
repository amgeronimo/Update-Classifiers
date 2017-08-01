function [viability, fileerr] = CheckFileViabilityP300(file,params)

fileerr = [];
viability = 1;

[a b c d] = load_bcidat(file,'-calibrated');
file;


%Do basic parameters exist?
try
    fs = c.SamplingRate.NumericValue;
    NumChans = size(a,2);
    NumSequences = c.NumberOfSequences.NumericValue;
    NumStimCodes = max(b.StimulusCode);
    OnTime = c.StimulusDuration.NumericValue/1000;
    OffTime = c.ISIMaxDuration.NumericValue/1000;
    SampleBlock = c.SampleBlockSize.NumericValue;
catch
    fileerr = [fileerr '~par_tsk'];
    viability = 0;
    return;
end


%These arent actually required, but files that have different values for these dont work together
if SampleBlock~=params.sbs
    fileerr = [fileerr '~par_sbs'];
    viability = 0;
end
if fs~=params.fs
    fileerr = [fileerr '~par_fs'];
    viability = 0;
end
if NumChans~=params.numch
    fileerr = [fileerr '~par_nc'];
    viability = 0;
end
if isfield(c,'ToBeCopied');
    tsc = 8;
else
    tsc = 18;
end
if NumStimCodes~=tsc
    fileerr = [fileerr '~par_nsc'];
    viability = 0;
end
if  OnTime ~= params.timeon
    fileerr = [fileerr '~par_ont'];
    viability = 0;
end
if  OffTime ~= params.timeoff
    fileerr = [fileerr '~par_oft'];
    viability = 0;
end

%Was it Copy Spelling?
if c.InterpretMode.NumericValue ~= 2
    fileerr = [fileerr '~fr'];
    viability = 0;
    return;
end


if isfield(c,'ToBeCopied') %Using audio speller
    T2S = c.ToBeCopied.NumericValue;
else %Using visual CB speller
    
    T2S = cell2mat(c.TextToSpell.Value);
    %Check for log file
    rstart = regexp(file, 'R+[0-9]+[(.\)]+');
    rend = regexp(file, '[(.\)]+dat');
    runnum = num2str(str2num(file(rstart+1:rend-1)));
    
    
    fid = fopen([file(1:rstart) runnum '+_mylogfile.txt'],'r');
    if fid == -1
        fileerr = [fileerr '~lg_mis'];
        viability = 0;
        return;
    else
        itt=textscan(fid,'%s','delimiter','\n');
        oldfile = sum(cellfun(@(x) ~isempty(strfind(x,'*** Beginning New Trial ***')), itt{:}));
        if oldfile>0
            fileerr = [fileerr '~lg_old'];
            viability = 0;
            return;
        else %Check if the file ends with two repetitions of 'Goal Text'
            GoalLind = find(cellfun(@(x) ~isempty(strfind(x,'Goal Text,')), itt{:}));
            ChosenLind = find(cellfun(@(x) ~isempty(strfind(x,'Selected Text')), itt{:}));
            trialkeep = ~ismember(GoalLind,ChosenLind+1);
            trialkeep(end)=0;
            GoalText = GoalLind(trialkeep);
            INITind = find(cellfun(@(x) ~isempty(strfind(x,'INIT_mSequence')), itt{:}));
            if length(INITind)==1 %see test001R78?
                fileerr = [fileerr '~lg_fin'];
                viability = 0;
                return;
            end
            if length(GoalText)~=length(regexp(T2S,'[A-Z1-9]'))
                fileerr = [fileerr '~lg_fin2'];
                viability = 0;
                return;
            end
            endlog = itt{1}(GoalLind(end-1:end));
            stidx = cellfun(@(x) regexp(x,'=\s'), endlog);
            endlog = cellfun(@(x) x(stidx+2:end),endlog,'UniformOutput',false);
            if sum(strcmp(endlog,T2S))~=2
                fileerr = [fileerr '~lg_fin3'];
                viability = 0;
                return;
            end
        end
        fclose(fid);
    end
end

%Is the data file complete?
remainder = mod(OnTime*fs,SampleBlock);
trimMax = round((c.PreRunDuration.NumericValue-.5)*fs);
StimulusCode = b.StimulusCode(trimMax:end);
temp=StimulusCode(StimulusCode>0);
temp=temp(1:OnTime*fs-remainder:end);
unique(temp);
for i = unique(temp)'
    NumTrials(i) = sum(temp==i)/NumSequences;
end


if length(temp)~=(NumSequences*NumStimCodes*mean(NumTrials)) || mean(NumTrials)~=NumTrials(1)
    fileerr = [fileerr '~dt'];
    viability = 0;
end


%If Data viable, give it a discriminibility score (or other metric to
%quickly rank for suitability in the classifier
if viability==1
    [bb,aa] = butter(3,[1 20]/(fs/2));
    a = filtfilt(bb,aa,a);
    
    rng = (round(.45*fs):round(1.2*fs));
    lrng = length(rng);
    targ = find(b.StimulusCode > 0 & b.StimulusType == 1);
    targ = targ(1:c.SampleBlockSize.NumericValue:end);
    targ = targ+rng(1);
    nontarg = find(b.StimulusCode>0 & b.StimulusType == 0);
    nontarg = nontarg(1:c.SampleBlockSize.NumericValue:end);
    nontarg = nontarg+rng(1);
    sdiff = zeros(lrng,4); kk=1;
    for i = [1 2 4 8]
        targ_d = a((repmat(targ,1,lrng)+repmat(0:lrng-1,length(targ),1))',i);
        targ_d = reshape(targ_d,lrng,length(targ_d)/lrng);
        %Remove outliers
        keeptrl = (sum(targ_d>75)+sum(targ_d<-75))==0;
        targ_d = targ_d(:,keeptrl);
        nontarg_d = a((repmat(nontarg,1,lrng)+repmat(0:lrng-1,length(nontarg),1))',i);
        nontarg_d = reshape(nontarg_d,lrng,length(nontarg_d)/lrng);
        %Remove outliers
        keeptrl = (sum(nontarg_d>75)+sum(nontarg_d<-75))==0;
        nontarg_d = nontarg_d(:,keeptrl);
        sdiff(:,kk) = (mean(targ_d,2)-mean(nontarg_d,2))./(std(targ_d,[],2).^2+std(nontarg_d,[],2).^2);
        kk=kk+1;
    end
    
    viability = mean(mean(abs(sdiff)));
    
    %     figure
    %     subplot(1,2,1); hold on; plot(rng/fs,nontarg_d,'Color',[1 .8 .8],'LineWidth',.5);
    %     plot(rng/fs,targ_d,'Color',[.8 .8 1],'LineWidth',.5);
    %     plot(rng/fs,mean(targ_d,2),'b','LineWidth',2);
    %     plot(rng/fs,mean(nontarg_d,2),'r','LineWidth',2);
    %     title(file(end-10:end-4))
    %     subplot(1,2,2);plot(sdiff); ylim([-.2 .2]);
    %     title(['v-score = ' num2str(viability)]);
end

end
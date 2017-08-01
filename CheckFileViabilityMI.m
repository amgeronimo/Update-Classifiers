function [viability, fileerr] = CheckFileViabilityMI(file,params)

fileerr = [];
viability = 1;

[a b c d] = load_bcidat(file,'-calibrated');


%Do basic parameters exist?
try
    fs = c.SamplingRate.NumericValue;
    NumChans = size(a,2);
    NumTarg = c.NumberTargets.NumericValue;
    SampleBlock = c.SampleBlockSize.NumericValue;
    NumTrials = c.NumberOfTrials.NumericValue;
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
if NumTarg~=2
    fileerr = [fileerr '~par_trg'];
    viability = 0;
end


%Is the data file complete?
TargCode = b.TargetCode;
trimMax = round((c.PreRunDuration.NumericValue-.5)*fs);
TargCode(1:trimMax) = 0;
temp=diff(TargCode);
temp=temp(temp~=0);
if length(temp)~=NumTrials
    fileerr = [fileerr '~dt'];
    viability = 0;
end

% 
% %If Data viable, give it a discriminibility score (or other metric to
% %quickly rank for suitability in the classifier
% if viability==1
%     [bb,aa] = butter(3,[1 20]/(fs/2));
%     a = filtfilt(bb,aa,a);
%     
%     rng = (round(.45*fs):round(1.2*fs));
%     lrng = length(rng);
%     targ = find(b.StimulusCode > 0 & b.StimulusType == 1);
%     targ = targ(1:c.SampleBlockSize.NumericValue:end);
%     targ = targ+rng(1);
%     nontarg = find(b.StimulusCode>0 & b.StimulusType == 0);
%     nontarg = nontarg(1:c.SampleBlockSize.NumericValue:end);
%     nontarg = nontarg+rng(1);
%     sdiff = zeros(lrng,4); kk=1;
%     for i = [1 2 4 8]
%         targ_d = a((repmat(targ,1,lrng)+repmat(0:lrng-1,length(targ),1))',i);
%         targ_d = reshape(targ_d,lrng,length(targ_d)/lrng);
%         %Remove outliers
%         keeptrl = (sum(targ_d>75)+sum(targ_d<-75))==0;
%         targ_d = targ_d(:,keeptrl);
%         nontarg_d = a((repmat(nontarg,1,lrng)+repmat(0:lrng-1,length(nontarg),1))',i);
%         nontarg_d = reshape(nontarg_d,lrng,length(nontarg_d)/lrng);
%         %Remove outliers
%         keeptrl = (sum(nontarg_d>75)+sum(nontarg_d<-75))==0;
%         nontarg_d = nontarg_d(:,keeptrl);
%         sdiff(:,kk) = (mean(targ_d,2)-mean(nontarg_d,2))./(std(targ_d,[],2).^2+std(nontarg_d,[],2).^2);
%         kk=kk+1;
%     end
%     
%     viability = mean(mean(abs(sdiff)));
%     
%     %     figure
%     %     subplot(1,2,1); hold on; plot(rng/fs,nontarg_d,'Color',[1 .8 .8],'LineWidth',.5);
%     %     plot(rng/fs,targ_d,'Color',[.8 .8 1],'LineWidth',.5);
%     %     plot(rng/fs,mean(targ_d,2),'b','LineWidth',2);
%     %     plot(rng/fs,mean(nontarg_d,2),'r','LineWidth',2);
%     %     title(file(end-10:end-4))
%     %     subplot(1,2,2);plot(sdiff); ylim([-.2 .2]);
%     %     title(['v-score = ' num2str(viability)]);
% end

end
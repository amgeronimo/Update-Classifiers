currdir = cd;


%Data recorded with the gUSBamp, CrossFace paradigm
cd('C:\Users\ageronimo\Google Drive\Hershey_2015\Study40647_2015\Classifiers\P300Classifier')
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD,ClassLabel, ChanLabel, SC, SR] = ...
    BuildClassifierFromTrainingData_v3('Andrew',{'010'},{{'01','02','03'}},1,1,0,[],[],[],[],[]);
%saved as gUSBLag.png


%Data recorded with the gNautilus, CrossFace paradigm
cd(currdir);
[nCM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    nTrainD,nClassLabel, nChanLabel, SC, SR] = ...
    BuildClassifier_teleBCI('WetTestAG',{'001'},{{'01','02','04'}},1,1,0,[],[],[],[],[]);
nChannelNames = {'Fp1','Fp2','F3','Fz','F4','T7','C3','Cz','C4','T8',...
        'P3','Pz','P4','PO7','PO8','Oz'};
%saved as gNautilusLag.png




%%  Tested varying blocktime (DryTestAG003)
%For this to work, temporarily change ReducT to 0 in 
%BuildClassifier_teleBCI.m (should normally be 3).  Al

%Blocktime of 25
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD1,ClassLabel1, ChanLabel1, SC, SR] = ...
    BuildClassifier_teleBCI('DryTestAG',{'003'},{{'01'}},1,1,0,[],[],[],[],[]);

%Blocktime of 12
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD2,ClassLabel2, ChanLabel2, SC, SR] = ...
    BuildClassifier_teleBCI('DryTestAG',{'003'},{{'11'}},1,1,0,[],[],[],[],[]);    
  
%Blocktime of 8
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD3,ClassLabel3, ChanLabel3, SC, SR] = ...
    BuildClassifier_teleBCI('DryTestAG',{'003'},{{'13'}},1,1,0,[],[],[],[],[]);
ChannelNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};

fs = 250;
Times = 1:fs;


    
figure
for i = 1:8
    subplot(2,4,i); hold on
    plot(Times/fs,mean(TrainD1(ClassLabel1(ChanLabel1==1)==1,(i-1)*fs+1:i*fs)));
    plot(Times/fs,mean(TrainD2(ClassLabel2(ChanLabel2==1)==1,(i-1)*fs+1:i*fs)));
    plot(Times/fs,mean(TrainD3(ClassLabel3(ChanLabel3==1)==1,(i-1)*fs+1:i*fs)));
end
legend('25','12','8')

%For Channels Fz, Cz,P4, PO8, Oz (1,2,5,7, 8 rows) peak latency for 
%R03,R11,R13 (cols)


pks = [.672 .44 .384; .672 .46 .384; .64 .404 .336; .66 .404 .356; .636 .4 .352];
diff(pks,[],2)
    
figure
barweb_AG(mean(pks),std(pks))

%%  Tested blocktime again (DryTestAG005)
%For this to work, temporarily change ReducT to 0 in 
%BuildClassifier_teleBCI.m (should normally be 3).  Al

%Blocktime of 25
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD1,ClassLabel1, ChanLabel1, SC, SR] = ...
    BuildClassifier_teleBCI('DryTestAG',{'005'},{{'01','03',}},1,1,0,[],[],[],[],[]);

%Blocktime of 12
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD2,ClassLabel2, ChanLabel2, SC, SR] = ...
    BuildClassifier_teleBCI('DryTestAG',{'005'},{{'04','06'}},1,1,0,[],[],[],[],[]);    
  
%Blocktime of 8
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD3,ClassLabel3, ChanLabel3, SC, SR] = ...
    BuildClassifier_teleBCI('DryTestAG',{'005'},{{'07','08'}},1,1,0,[],[],[],[],[]);
ChannelNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};

fs = 250;
Times = 1:fs;


    
figure
for i = 1:8
    subplot(2,4,i); hold on
    plot(Times/fs,mean(TrainD1(ClassLabel1(ChanLabel1==1)==1,(i-1)*fs+1:i*fs)));
    plot(Times/fs,mean(TrainD2(ClassLabel2(ChanLabel2==1)==1,(i-1)*fs+1:i*fs)));
    plot(Times/fs,mean(TrainD3(ClassLabel3(ChanLabel3==1)==1,(i-1)*fs+1:i*fs)));
end
legend('25','12','8')

%For Channels Fz, Cz,P3, Pz, P4, PO8 (1,2,3,4,5,7) peak latency for 
%the three block sizes( cols)


pks = [.66 .44 .368;
    .664 .444 .368;
    .68 .46 .392;
    .66 .404 .356;
    .684 .456 .408;
    .716 .504 .428;
    .716 .504 .436];
actpk = mean(diff(pks,[],2))
theorpk = [8/250*4-25/250*4 8/250*4-12/250*4]
    
figure
barweb_AG(mean(pks),std(pks))

%%  Tested blocktime again (WetTestAG005)
%For this to work, temporarily change ReducT to 0 in 
%BuildClassifier_teleBCI.m (should normally be 3).  Al

%Blocktime of 25
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD1,ClassLabel1, ChanLabel1, SC, SR] = ...
    BuildClassifier_teleBCI('WetTestAG',{'005'},{{'01','02',}},1,1,0,[],[],[],[],[]);

%Blocktime of 12
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD2,ClassLabel2, ChanLabel2, SC, SR] = ...
    BuildClassifier_teleBCI('WetTestAG',{'005'},{{'03','04'}},1,1,0,[],[],[],[],[]);    
  
%Blocktime of 8
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD3,ClassLabel3, ChanLabel3, SC, SR] = ...
    BuildClassifier_teleBCI('WetTestAG',{'005'},{{'05','06'}},1,1,0,[],[],[],[],[]);
nChannelNames = {'Fp1','Fp2','F3','Fz','F4','T7','C3','Cz','C4','T8',...
        'P3','Pz','P4','PO7','PO8','Oz'};

fs = 250;
Times = 1:fs;


    
figure
for i = 1:16
    subplot(4,4,i); hold on
    plot(Times/fs,mean(TrainD1(ClassLabel1(ChanLabel1==1)==1,(i-1)*fs+1:i*fs)));
    plot(Times/fs,mean(TrainD2(ClassLabel2(ChanLabel2==1)==1,(i-1)*fs+1:i*fs)));
    plot(Times/fs,mean(TrainD3(ClassLabel3(ChanLabel3==1)==1,(i-1)*fs+1:i*fs)));
end
legend('25','12','8')

%For Channels 1,2,3,4,5,6,7,8,13,14,15 peak latency for 
%the three block sizes( cols)


pks = [.656 .448 .392;
    .66 .46 .388;
    .664 .452 .392;
    .664 .452 .392;
    .664 .452 .392;
    .652 .452 .396;
    .668 .452 .392;
    .668 .46 .392;
    .632 .416 .36;
    .632 .42 .356;
    .644 .436 .38];
actpk = mean(diff(pks,[],2))
theorpk = [8/250*4-25/250*4 8/250*4-12/250*4]
    
figure
barweb_AG(mean(pks),std(pks))
%% Check out the sourcetime and stimulus time parameters
%25 block
[a1,b1,c1,d1] = load_bcidat(['C:\Users\ageronimo\Documents\BCI2000_5300\'...
    'data\test001\testS001R90.dat'],'-calibrated');
%8 block
[a2,b2,c2,d2] = load_bcidat(['C:\Users\ageronimo\Documents\BCI2000_5300\'...
    'data\test001\testS001R91.dat'],'-calibrated');
%1 block
[a3,b3,c3,d3] = load_bcidat(['C:\Users\ageronimo\Documents\BCI2000_5300\'...
    'data\test001\testS001R92.dat'],'-calibrated');



figure
xx(1) = subplot(1,3,1); plot(b1.SourceTime); hold on; plot(b1.StimulusTime);
legend('SourceTime','StimulusTime'); title('Blocksize 25');
xx(2) = subplot(1,3,2); plot(b2.SourceTime); hold on; plot(b2.StimulusTime);
title('Blocksize 8');
xx(3) = subplot(1,3,3); plot(b3.SourceTime); hold on; plot(b3.StimulusTime);
title('Blocksize 1');
linkaxes(xx,'x');




%% Do the same think but with the signalgenerator
%25 block
[a1,b1,c1,d1] = load_bcidat(['C:\Users\ageronimo\Documents\BCI2000_5300\'...
    'data\test001\testS001R94.dat'],'-calibrated');
%8 block
[a2,b2,c2,d2] = load_bcidat(['C:\Users\ageronimo\Documents\BCI2000_5300\'...
    'data\test001\testS001R95.dat'],'-calibrated');
%1 block
[a3,b3,c3,d3] = load_bcidat(['C:\Users\ageronimo\Documents\BCI2000_5300\'...
    'data\test001\testS001R96.dat'],'-calibrated');



figure
xx(1) = subplot(1,3,1); plot(b1.SourceTime); hold on; plot(b1.StimulusTime);
legend('SourceTime','StimulusTime'); title('Blocksize 25');
xx(2) = subplot(1,3,2); plot(b2.SourceTime); hold on; plot(b2.StimulusTime);
title('Blocksize 8');
xx(3) = subplot(1,3,3); plot(b3.SourceTime); hold on; plot(b3.StimulusTime);
title('Blocksize 1');
linkaxes(xx,'x');

%For these to run, need to set channels 1:16
%Blocktime of 25
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD1,ClassLabel1, ChanLabel1, SC, SR] = ...
    BuildClassifier_teleBCI('test',{'001'},{{'94'}},1,1,0,[],[],[],[],[]);
%Blocktime of 8
[CM, allsums, targetstim, WORD_Acc, skiptest, TTSpell, ChosenLetterF,...
    TrainD3,ClassLabel3, ChanLabel3, SC, SR] = ...
    BuildClassifier_teleBCI('test',{'001'},{{'95'}},1,1,0,[],[],[],[],[]);
%Blocktime of 1 (did not finish, cannot run)


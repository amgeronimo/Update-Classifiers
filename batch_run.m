%Do batch runs of data for P300
%  [CM,~,~,WORDD, skiptest] = BuildClassifierFromTrainingData('','',{''},1,1,0,[],[],[],[],[]);
% Filename = {'Andrew
% 
% 
% for ff = 1:size(Filename,1)
% filename = 'AndrewP3_14';
% filesess = '001';
% filerun = '07';

% clear all
close all

%Copy for non-batch runs
%[CM,~,~,WORDD, skiptest] = BuildClassifierFromTrainingData('','',{''},1,1,0,[],[],[],[],[]);

% filename = 'S09_CrossCB_22_3EOG_2amp';
% filename = 'C02';]
% filesess = {'001','002'};
% filerun = {{'11','12','13','14'},{'01','02'}};

filename = 'WetTestAG';
filesess = {'001'};
filerun = {{'01','02','04'}};

%All of my tests show that the SWLDA and regularized LDA (lasso) work fine
Rerefs = [1];% 1 2 3];
CSPs = [0];
Arts = [0];
Classifs = [1]; % 5( RegLDA cannot be run on Matlab2012

totalcombos = length(Rerefs)*length(CSPs)*length(Arts)*length(Classifs);

wtbrr = waitbar(0,'Progress');
totalind = 1;

clear CM WORDD T_WORDD T_WORDD2 T_WORDD3 T_WORDD4 T_WORDD5
r_ind = 1;
for RV = Rerefs
    d_ind = 1;
    for DCSP = CSPs
        a_ind = 1;
        for ART = Arts
            c_ind = 1;
            for Cfier = Classifs
                %Build the classifier
                [~,~,~,WORDD(r_ind,d_ind,a_ind,c_ind), skiptest] = ...
                    BuildClassifier_teleBCI(filename,...
                    filesess,filerun,1,Cfier,1,RV,ART,DCSP,[],[]);
                
                if skiptest == 0
                    testsess = {'001'};
                    testrun = {{'09'}};
                    [CM{r_ind,d_ind,c_ind},~,~,tWORDD] = ...
                        BuildClassifier_teleBCI(filename,...
                        filesess(end),testrun,0,Cfier,1,RV,ART,DCSP,testsess,filerun);
                    if isnumeric(tWORDD)
                        T_WORDD(r_ind,d_ind,a_ind,c_ind) = tWORDD;
                    else
                        T_WORDD{r_ind,d_ind,a_ind,c_ind} = tWORDD;
                    end
%                     testsess = {'002'};
%                     testrun = {{'04'}};
%                     [CM2{r_ind,d_ind,c_ind},~,~,tWORDD] = ...
%                         BuildClassifierFromTrainingData_v3(filename,...
%                         filesess(end),testrun,0,Cfier,1,RV,ART,DCSP,testsess,filerun);
%                     if isnumeric(tWORDD)
%                         T_WORDD2(r_ind,d_ind,a_ind,c_ind) = tWORDD;
%                     else
%                         T_WORDD2{r_ind,d_ind,a_ind,c_ind} = tWORDD;
%                     end
%                     
%                     testrun = '07';
%                     [CM3{r_ind,d_ind,c_ind},~,~,tWORDD] = ...
%                         BuildClassifierFromTrainingData_v3(filename,...
%                         filesess,{testrun},0,Cfier,1,RV,ART,DCSP,filesess,filerun);
%                     if isnumeric(tWORDD)
%                         T_WORDD3(r_ind,d_ind,a_ind,c_ind) = tWORDD;
%                     else
%                         T_WORDD3{r_ind,d_ind,a_ind,c_ind} = tWORDD;
%                     end
%                     testrun = '06';
%                     [CM,~,~,T_WORDD4(r_ind,d_ind,a_ind,c_ind)] = ...
%                         BuildClassifierFromTrainingData(filename,...
%                         filesess,{testrun},0,Cfier,1,RV,ART,DCSP,filesess,filerun);
                    
%                     testrun = '06';
%                     [CM,~,~,T_WORDD5(r_ind,d_ind,a_ind,c_ind)] = ...
%                         BuildClassifierFromTrainingData(filename,...
%                         filesess,{testrun},0,Cfier,1,RV,ART,DCSP,filesess,filerun);
                    
                    
                end
                waitbar(totalind/totalcombos,wtbrr);
                totalind = totalind+1;
                c_ind = c_ind+1;
            end
            a_ind = a_ind+1;
        end
        d_ind = d_ind+1;
    end
    r_ind = r_ind+1;
end
           
close(wtbrr)
    

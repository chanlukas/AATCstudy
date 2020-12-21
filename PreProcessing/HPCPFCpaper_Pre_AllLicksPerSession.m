clear all

%% Load Data from HPConly Experiment
cd('T:\jan\Collabo Data')
load AATCstrucPFConly.mat
load AATCpfcAllBehavior.mat

for i=1:length(AATCpfcAllBehavior)  
BehaviorAnimals(i)=AATCpfcAllBehavior(i).Animal;
BehaviorSession(i)=AATCpfcAllBehavior(i).Session;
end

counter=1;
for i=1:length(AATC_PFC)
    i
 
LearnedCounter(counter)=AATC_PFC(i).Learned(1,1)

idx= find(BehaviorAnimals==(cell2mat(AATC_PFC(i).Animal)-7)&BehaviorSession==cell2mat(AATC_PFC(i).TrainingDay))
Session=AATCpfcAllBehavior(idx).Events;
if isempty(Session{1,1})~=1
pre=1*1000-1;
post=4*1000;
bin=100;
Fig=0;
[RewardMatrix,unRewardMatrix]=HPCpaper_Pre_BehaviorExample(Session,pre,post,bin,Fig);
AvgLickCounter(counter,:,1)=mean(RewardMatrix);
AvgLickCounter(counter,:,2)=mean(unRewardMatrix);
else
AvgLickCounter(counter,1:5000,1)=NaN;
AvgLickCounter(counter,1:5000,2)=NaN;   
end
counter=counter+1;

end

clearvars -except  LearnedCounter counter AvgLickCounter
%% Load Data from HPConly Experiment
cd('T:\jan\Collabo Data')
load AATCstrucHPConly.mat
load AATChpcAllBehavior.mat

for i=1:length(AATChpcAllBehavior)  
BehaviorAnimals(i)=AATChpcAllBehavior(i).Animal;
BehaviorSession(i)=AATChpcAllBehavior(i).Session;
end

counter=1;
for i=1:length(AATC_HPC)
    i


idx= find(BehaviorAnimals==cell2mat(AATC_HPC(i).Animal)&BehaviorSession==cell2mat(AATC_HPC(i).TrainingDay))
Session=AATChpcAllBehavior(idx).Events;
pre=1*1000-1;
post=4*1000;
bin=100;
Fig=0;
[RewardMatrix,unRewardMatrix]=HPCpaper_Pre_BehaviorExample(Session,pre,post,bin,Fig);
AvgLickCounter(counter,:,1)=mean(RewardMatrix);
AvgLickCounter(counter,:,2)=mean(unRewardMatrix);
LearnedCounter(counter)=AATC_HPC(i).Learned(1,1)


counter=counter+1;

end

clearvars -except  LearnedCounter counter AvgLickCounter

load  AATCstrucHPCPFC.mat
load  AATChpcpfcAllBehavior.mat

%%Exp2

for i=1:length(AATChpcpfcAllBehavior)  
BehaviorAnimals(i)=AATChpcpfcAllBehavior(i).Animal;
BehaviorSession(i)=AATChpcpfcAllBehavior(i).Session;
end

for i=1:length(AATC_HPC_PFC)


idx= find(BehaviorAnimals==(cell2mat(AATC_HPC_PFC(i).Animal)-13)&BehaviorSession==cell2mat(AATC_HPC_PFC(i).TrainingDay))
Session=AATChpcpfcAllBehavior(idx).Events;
pre=1*1000-1;
post=4*1000;
bin=100;
Fig=0;
[RewardMatrix,unRewardMatrix]=HPCpaper_Pre_BehaviorExample(Session,pre,post,bin,Fig);
AvgLickCounter(counter,:,1)=mean(RewardMatrix);
AvgLickCounter(counter,:,2)=mean(unRewardMatrix);

LearnedCounter(counter)=AATC_HPC_PFC(i).Learned(1,1)

counter=counter+1;
end

cd('T:\jan\Collabo Data\PFCpaperPreProcessed')

save HPCPFC_AllLicks AvgLickCounter LearnedCounter



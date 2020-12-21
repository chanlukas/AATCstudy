clear all

%% Load Data from HPConly Experiment
cd('T:\jan\Collabo Data')
load AATCstrucHPConly.mat
load AATChpcAllBehavior.mat

for i=1:length(AATChpcAllBehavior)  
BehaviorAnimals(i)=AATChpcAllBehavior(i).Animal;
BehaviorSession(i)=AATChpcAllBehavior(i).Session;
end

cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
load('LickEvokedIndx.mat')

counter=1;
Cellcounter=1;

counter=1;
for i=1:length(AATC_HPC)
    i
   SpikeTS=AATC_HPC(i).Spikes{1,1};
   ClusterID=AATC_HPC(i).Spikes{2,1};
   GoodClusters=unique(ClusterID);
   
   EventTypeAll=AATC_HPC(i).Events{2,1};
   EventTSAll=AATC_HPC(i).Events{1,1};
   
   EventType=EventTypeAll(ismember(EventTypeAll,[1,2]));
   EventTS=EventTSAll(ismember(EventTypeAll,[1,2]));
   
   RecordingLength=AATC_HPC(i).RecordingLength;
   
   Pre=1;
   Post=4;
   Bin=1;
   
[PsthSua,EventMatrix,EventType]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,Pre,Post,Bin);


idx1=find(LickUpPFC((Cellcounter:Cellcounter+size(PsthSua,2)-1))==0&LickDownPFC((Cellcounter:Cellcounter+size(PsthSua,2)-1))==0)
data=squeeze(sum(EventMatrix(:,90001:120000,idx1),2));
base=squeeze(sum(EventMatrix(:,1:15000,idx1),2));
Rwd=squeeze(sum(EventMatrix(:,120001:150000,idx1),2));

Cellcounter=Cellcounter+size(PsthSua,2);

if length(find(EventType==1))>length(find(EventType==2))
  idx=[randsample(1:length(find(EventType==1)),length(find(EventType==2))),...
       length(find(EventType==1))+1:length(EventType)];
else
      idx=[1:length(find(EventType==1)),...
          randsample(length(find(EventType==1))+1:length(EventType),length(find(EventType==1)))];
end

Y = num2cell(num2str(EventType(idx)));

MdlTrace = fitcsvm(data(idx,:),Y);
MdlBase = fitcsvm(base(idx,:),Y);
MdlRwd = fitcsvm(Rwd(idx,:),Y);


McRTrace = kfoldLoss(crossval(MdlTrace, 'KFold', 20))

McRBase = kfoldLoss(crossval(MdlBase, 'KFold', 20))

McRRwd = kfoldLoss(crossval(MdlRwd, 'KFold', 20))


allMcRT(counter,:)=McRTrace;
allMcRb(counter,:)=McRBase;
allMcRR(counter,:)=McRRwd;

LearnedCounter(counter)=AATC_HPC(i).Learned(1,1)

% idx= find(BehaviorAnimals==cell2mat(AATC_HPC(i).Animal)&BehaviorSession==cell2mat(AATC_HPC(i).TrainingDay))
% Session=AATChpcAllBehavior(idx).Events;
% pre=1*1000-1;
% post=4*1000;
% bin=100;
% Fig=0;
% [RewardMatrix,unRewardMatrix]=HPCpaper_Pre_BehaviorExample(Session,pre,post,bin,Fig);
% AvgLickCounter(counter,:,1)=mean(RewardMatrix);
% AvgLickCounter(counter,:,2)=mean(unRewardMatrix);


counter=counter+1;

end

clearvars -except allMcRb allMcRR allMcRT LearnedCounter counter Cellcounter AvgLickCounter LickUpPFC LickDownPFC
cd('T:\jan\Collabo Data')

load  AATCstrucHPCPFC.mat
load  AATChpcpfcAllBehavior.mat

%%Exp2

for i=1:length(AATChpcpfcAllBehavior)  
BehaviorAnimals(i)=AATChpcpfcAllBehavior(i).Animal;
BehaviorSession(i)=AATChpcpfcAllBehavior(i).Session;
end

for i=1:length(AATC_HPC_PFC)
    i
   SpikeTS=AATC_HPC_PFC(i).Spikes{1,1};
   ClusterID=AATC_HPC_PFC(i).Spikes{2,1};
   ShankS=cell2mat(AATC_HPC_PFC(i).Area);
   GoodClusters=unique(ClusterID);
   GoodClusters=GoodClusters(ShankS==1);
   
   EventTypeAll=AATC_HPC_PFC(i).Events{2,1};
   EventTSAll=AATC_HPC_PFC(i).Events{1,1};
   
   EventType=EventTypeAll(ismember(EventTypeAll,[1,2]));
   EventTS=EventTSAll(ismember(EventTypeAll,[1,2]));
   
   RecordingLength=AATC_HPC_PFC(i).RecordingLength;
   
   Pre=1;
   Post=5;
   Bin=1;
   
   [PsthSua,EventMatrix,EventType]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,Pre,Post,Bin);

   
idx1=find(LickUpPFC((Cellcounter:Cellcounter+size(PsthSua,2)-1))==0&LickDownPFC((Cellcounter:Cellcounter+size(PsthSua,2)-1))==0)

data=squeeze(sum(EventMatrix(:,90001:120000,idx1),2));
base=squeeze(sum(EventMatrix(:,1:15000,idx1),2));
Rwd=squeeze(sum(EventMatrix(:,120001:150000,idx1),2));

Cellcounter=Cellcounter+size(PsthSua,2);


if length(find(EventType==1))>length(find(EventType==2))
  idx=[randsample(1:length(find(EventType==1)),length(find(EventType==2))),...
       length(find(EventType==1))+1:length(EventType)];
else
      idx=[1:length(find(EventType==1)),...
          randsample(length(find(EventType==1))+1:length(EventType),length(find(EventType==1)))];
end

Y = num2cell(num2str(EventType(idx)));

MdlTrace = fitcsvm(data(idx,:),Y);
MdlBase = fitcsvm(base(idx,:),Y);
MdlRwd = fitcsvm(Rwd(idx,:),Y);


McRTrace = kfoldLoss(crossval(MdlTrace, 'KFold', 20))

McRBase = kfoldLoss(crossval(MdlBase, 'KFold', 20))

McRRwd = kfoldLoss(crossval(MdlRwd, 'KFold', 20))


allMcRT(counter,:)=McRTrace;
allMcRb(counter,:)=McRBase;
allMcRR(counter,:)=McRRwd;

LearnedCounter(counter)=AATC_HPC_PFC(i).Learned(1,1)

% idx= find(BehaviorAnimals==(cell2mat(AATC_HPC_PFC(i).Animal)-13)&BehaviorSession==cell2mat(AATC_HPC_PFC(i).TrainingDay))
% Session=AATChpcpfcAllBehavior(idx).Events;
% pre=1*1000-1;
% post=4*1000;
% bin=100;
% Fig=0;
% [RewardMatrix,unRewardMatrix]=HPCpaper_Pre_BehaviorExample(Session,pre,post,bin,Fig);
% AvgLickCounter(counter,:,1)=mean(RewardMatrix);
% AvgLickCounter(counter,:,2)=mean(unRewardMatrix);
% 
counter=counter+1;
end

cd('T:\jan\Collabo Data\HPCpaperPreProcessed')

save SVMLickCorrected allMcRb allMcRR allMcRT LearnedCounter AvgLickCounter


n1{1}=1-allMcRb(LearnedCounter==1)
n1{2}=1-allMcRT(LearnedCounter==1)
names={'Baseline Learned','Trace Learned'};
ttl={'SVM Classifier'};
colors=[0 1 0;0 0 0;];
ylabel={'Percent Correct'}
style=([1 1 1;1 1 1]);
figure()
BeehivePlot(n1,names,ylabel,ttl,colors,style)
hold on
plot([.5 2.5],[0.5 0.5],':r')
set(gca,'FontSize',25);
set(gca,'LineWidth',5);



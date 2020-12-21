clear all

%% Load Data from HPConly Experiment
cd('T:\jan\Collabo Data')
load AATCstrucPFConly.mat
load AATCpfcAllBehavior.mat

for i=1:length(AATCpfcAllBehavior)  
BehaviorAnimals(i)=AATCpfcAllBehavior(i).Animal;
BehaviorSession(i)=AATCpfcAllBehavior(i).Session;
end

cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
load('LickEvokedIndx.mat')

counter=1;
Cellcounter=1;

for i=1:length(AATC_PFC)
    i
   SpikeTS=AATC_PFC(i).Spikes{1,1};
   ClusterID=AATC_PFC(i).Spikes{2,1};
   GoodClusters=unique(ClusterID);
   
   EventTypeAll=AATC_PFC(i).Events{2,1};
   EventTSAll=AATC_PFC(i).Events{1,1};
   
   EventType=EventTypeAll(ismember(EventTypeAll,[1,2]));
   EventTS=EventTSAll(ismember(EventTypeAll,[1,2]));
   
   RecordingLength=AATC_PFC(i).RecordingLength;
   
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

LearnedCounter(counter)=AATC_PFC(i).Learned(1,1)

% idx= find(BehaviorAnimals==(cell2mat(AATC_PFC(i).Animal)-7)&BehaviorSession==cell2mat(AATC_PFC(i).TrainingDay))
% Session=AATCpfcAllBehavior(idx).Events;
% if isempty(Session{1,1})~=1
% pre=1*1000-1;
% post=4*1000;
% bin=100;
% Fig=0;
% [RewardMatrix,unRewardMatrix]=HPCpaper_Pre_BehaviorExample(Session,pre,post,bin,Fig);
% AvgLickCounter(counter,:,1)=mean(RewardMatrix);
% AvgLickCounter(counter,:,2)=mean(unRewardMatrix);
% else
% AvgLickCounter(counter,1:5000,1)=NaN;
% AvgLickCounter(counter,1:5000,2)=NaN;   
% end
counter=counter+1;
clearvars PsthSua EventMatrix EventType
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
   GoodClusters=GoodClusters(ShankS==2);
   
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

   Cellcounter=Cellcounter+size(PsthSua,2);
   
data=squeeze(sum(EventMatrix(:,90001:120000,idx1),2));
base=squeeze(sum(EventMatrix(:,1:15000,idx1),2));
Rwd=squeeze(sum(EventMatrix(:,120001:150000,idx1),2));


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

counter=counter+1;
end

cd('T:\jan\Collabo Data\PFCpaperPreProcessed')

save SVMLickCorrected allMcRb allMcRR allMcRT LearnedCounter 


n1{1}=1-allMcRbPFC
n1{2}=1-allMcRTPFC
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



clear all
%% Load Data from HPConly Experiment
cd('T:\jan\Collabo Data')
load AATCR2strucHPCPFC.mat
load AATCR2hpcpfcAllBehavior.mat 


AATCR2_Sua_Psth=[];
wfs=[];
includeWF=[];
SpikeV2P=[];
spikeWidth=[];
Depth=[];
Shank=[];
FRperiodsCounter=[];
SessionCounter=1;
CellCounter=0;

for i=1:length(AATCR2Behavior)  
BehaviorAnimals(i)=AATCR2Behavior(i).Animal;
BehaviorSession(i)=AATCR2Behavior(i).Session;
end
for i=1:length(AATCR2)
    i
   SpikeTS=AATCR2(i).Spikes{1,1};
   ClusterID=AATCR2(i).Spikes{2,1};
   ShankS=cell2mat(AATCR2(i).Area);
   GoodClusters=unique(ClusterID);
%    GoodClusters=GoodClusters(ShankS==1);
   
   EventTypeAll=AATCR2(i).Events{1,1}{2,1};
   EventTSAll=AATCR2(i).Events{1,1}{1,1};
   
   EventType=EventTypeAll(ismember(EventTypeAll,[1,2,6]));
   EventTS=EventTSAll(ismember(EventTypeAll,[1,2,6]));
   
   RecordingLength=AATCR2(i).RecordingLength;
   
   Pre=1;
   Post=5;
   Bin=1;
   
   [PsthSua,EventMatrix,EventType]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,Pre,Post,Bin);

%% Anova  
Cells=1:length(GoodClusters);
        for ii=1:length(Cells)
        Baseline=squeeze(mean(EventMatrix(:,1:30000,ii),2));
        Evoked=squeeze(mean(EventMatrix(:,30000:45000,ii),2));
        Sustained=squeeze(mean(EventMatrix(:,45000:90000,ii),2));
        Trace=squeeze(mean(EventMatrix(:,90000:120000,ii),2));
        Reward=squeeze(mean(EventMatrix(:,120000:150000,ii),2));
        Data=cat(1,Baseline,Evoked,Sustained,Trace,Reward);
        TimeGroup=repelem(1:5,size(EventMatrix,1));
        if size(EventType,2)~=1
            EventType=EventType',
        end
        EventTypeGroup=repmat(EventType,5,1);
        [p] = anovan(Data,{TimeGroup',EventTypeGroup},'model','interaction','display','off');
        Anovap(ii,:)=p;
        if min(Anovap)<0.05
        TaskMod(ii)=1;
        else
        TaskMod(ii)=0;
        end

        end
   
%% Firing Rate
bins=10;
sets=floor(AATCR2(i).RecordingLength/10); 
for ii=1:length(GoodClusters)   
TScell(1:RecordingLength)=0;
TScell(SpikeTS(find(ClusterID==GoodClusters(ii))))=1;
Periods=reshape(TScell(1:sets*(bins)),(bins),sets);
FRperiods(ii,:)=(sum(Periods,2)/(sets/30000))';
FRmean(ii)=mean(FRperiods(ii,:));
FRstd(ii)=std(FRperiods(ii,:));
end
% Licks

idx= find(BehaviorAnimals==cell2mat(AATCR2(i).Animal)&BehaviorSession==cell2mat(AATCR2(i).TrainingDay))
Session=AATCR2Behavior(idx).Events;
pre=1*1000-1;
post=5*1000;
bin=100;
Fig=0;
[RewardMatrix,RewardMatrix2,unRewardMatrix]=HPCpaper_Pre_BehaviorExampleR2(Session,pre,post,bin);
[Learned,LearnedP]=ttest2(mean(RewardMatrix(:,pre+2000:pre+2975),2),mean(unRewardMatrix(:,pre+2000:pre+2975),2));   
[Learned2,LearnedP2]=ttest2(mean(RewardMatrix2(:,pre+2000:pre+2975),2),mean(unRewardMatrix(:,pre+2000:pre+2975),2));   


LickRateDif=sum(sum(RewardMatrix(:,pre:pre+2975)))-sum(sum(unRewardMatrix(:,pre:pre+2975)));
LickRateR=mean(sum(RewardMatrix(:,pre:pre+2975)));
LickRateU=mean(sum(unRewardMatrix(:,pre:pre+2975)));
%% Counters
% Main Psths
AATCR2_Sua_Psth=cat(2,AATCR2_Sua_Psth,PsthSua);
    
TrgDayCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATCR2(i).TrainingDay);
AnimalCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATCR2(i).Animal);
SessionCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATCR2(i).Session);
LearnedCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=Learned;
LearnedPCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LearnedP;
Learned2Counter(CellCounter+1:CellCounter+size(GoodClusters,1))=Learned2;
Learned2PCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LearnedP2;
% CellPerSesCounter(CellCounter+1:CellCounter+size(PsthSua,2))=1:size(PsthSua,2);

LickRateDifCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LickRateDif;
LickRateRCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LickRateR;
LickRateUCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LickRateU;
AnovaCounter(CellCounter+1:CellCounter+size(PsthSua,2),:)=Anovap;
TaskModCounter(CellCounter+1:CellCounter+size(PsthSua,2))=TaskMod;
FRperiodsCounter=cat(1,FRperiodsCounter,FRperiods);

%%Waveform Extraction
wfsS=cell2mat(AATCR2(i).Waveforms.wfs);
ShankS=cell2mat(AATCR2(i).Area);
includeWFS=cell2mat(AATCR2(i).Waveforms.include);
SpikeV2PS=cell2mat(AATCR2(i).Waveforms.SpikeV2P);
SpikeWidthS=cell2mat(AATCR2(i).Waveforms.SpikeWidth);
DepthS=cell2mat(AATCR2(i).Waveforms.Depth);

%%Waveform Counters
wfs=cat(1,wfs,wfsS(ShankS==1,:));
includeWF=cat(2,includeWF,includeWFS);
SpikeV2P=cat(1,SpikeV2P,SpikeV2PS);
spikeWidth=cat(2,spikeWidth,SpikeWidthS);
Depth=cat(2,Depth,DepthS) ; 
Shank=cat(2,Shank,ShankS)

CellCounter=CellCounter+size(GoodClusters,1);
clearvars FRperiods PsthSua EventMatrix EventType wfsS includeWFS LickRateDif LickRateR LickRateU SpikeV2PS DepthS Anovap TaskMod ShankS

end
%% Waveform Analysis
GMMdata=horzcat(SpikeV2P(includeWF==1),spikeWidth(includeWF==1)',mean(FRperiodsCounter(includeWF==1,:),2))
GMMdata(GMMdata==0)=1;

GMModel = fitgmdist(GMMdata,2,'SharedCovariance',true)
[neuronTypeGMM,nlogL,P]  = cluster(GMModel,GMMdata);
neuronTypeGMM(P(:,1)>0.05&P(:,2)>0.05)=3;
neuronTypeGmm=zeros(1,length(AnimalCounter));
neuronTypeGmm(includeWF==1)=neuronTypeGMM;
figure()
scatter3(SpikeV2P(includeWF==1),spikeWidth(includeWF==1)',mean(FRperiodsCounter(includeWF==1,:),2),'g')
hold on
scatter3(SpikeV2P(neuronTypeGmm==1),spikeWidth(neuronTypeGmm==1)',mean(FRperiodsCounter(neuronTypeGmm==1,:),2),'r')
hold on
scatter3(SpikeV2P(neuronTypeGmm==2),spikeWidth(neuronTypeGmm==2)',mean(FRperiodsCounter(neuronTypeGmm==2,:),2),'b')

% Manual Threshold
neuronType(1:length(spikeWidth))=0;
neuronType(intersect(find(includeWF==1),find(SpikeV2P<=18)))=1;
neuronType(intersect(find(includeWF==1),find(SpikeV2P>18)))=2;

%% Save
cd('T:\jan\Collabo Data\HPCpaperPreProcessed')
save AATCR2_Sua_Psth_1ms AATCR2_Sua_Psth TrgDayCounter AnimalCounter ...
                    CellCounter SessionCounter...
                    wfs SpikeV2P spikeWidth neuronTypeGmm Depth includeWF ...
                    AnovaCounter TaskModCounter ...
                    FRperiodsCounter neuronType Shank ...
                    LickRateDifCounter LearnedCounter LearnedPCounter Learned2Counter Learned2PCounter LickRateRCounter LickRateUCounter


clear all

%% Load Data from HPConly Experiment
cd('T:\jan\Collabo Data')
load AATCstrucPFConly.mat
load AATCpfcAllBehavior.mat
AATC_Sua_Psth=[];
wfs=[];
includeWF=[];
SpikeV2P=[];
spikeWidth=[];
Depth=[];
FRperiodsCounter=[];
SessionCounter=1;
CellCounter=0;
%%Exp1
for i=1:length(AATCpfcAllBehavior)  
BehaviorAnimals(i)=AATCpfcAllBehavior(i).Animal;
BehaviorSession(i)=AATCpfcAllBehavior(i).Session;
end
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
 
   PFCCells=find(cell2mat(AATC_PFC(i).Area)==2);
   for ii=1:length(PFCCells)
   EvokedSig(ii)=ranksum(squeeze(mean(EventMatrix(EventType==1,Pre*30000:(Pre+1.35)*30000,PFCCells(ii)),2)),...
   squeeze(mean(EventMatrix(EventType==2,Pre*30000:(Pre+1.35)*30000,PFCCells(ii)),2))) ;

   TraceSig(ii)=ranksum(squeeze(mean(EventMatrix(EventType==1,(Pre+2)*30000:(Pre+3)*30000,PFCCells(ii)),2)),...
   squeeze(mean(EventMatrix(EventType==2,(Pre+2)*30000:(Pre+3)*30000,PFCCells(ii)),2))) ;
   end
% %% Anova  
        for ii=1:length(PFCCells)
        Baseline=squeeze(mean(EventMatrix(:,1:Pre*30000,PFCCells(ii)),2));
        Evoked=squeeze(mean(EventMatrix(:,Pre*30000:Pre*30000+15000,PFCCells(ii)),2));
        Sustained=squeeze(mean(EventMatrix(:,Pre*30000+15000:Pre*30000+60000,PFCCells(ii)),2));
        Trace=squeeze(mean(EventMatrix(:,Pre*30000+60000:Pre*30000+90000,PFCCells(ii)),2));
        Reward=squeeze(mean(EventMatrix(:,Pre*30000+90000:Pre*30000+120000,PFCCells(ii)),2));
        Data=cat(1,Baseline,Evoked,Sustained,Trace,Reward);
        TimeGroup=repelem(1:5,size(EventMatrix,1));
        EventTypeGroup=repmat(EventType,5,1);
        [p,tbl,stats] = anovan(Data,{TimeGroup',EventTypeGroup},'model','interaction','display','off');
        Anovap(ii,:)=p;
        if min(Anovap)<0.05
        TaskMod(ii)=1;
        else
        TaskMod(ii)=0;
        end
        end
   
%% Firing Rate
bins=10;
sets=floor(AATC_PFC(i).RecordingLength/10); 
for ii=1:length(GoodClusters)   
TScell(1:RecordingLength)=0;
TScell(SpikeTS(find(ClusterID==GoodClusters(ii))))=1;
Periods=reshape(TScell(1:sets*(bins)),(bins),sets);
FRperiods(ii,:)=(sum(Periods,2)/(sets/30000))';
FRmean(ii)=mean(FRperiods(ii,:));
FRstd(ii)=std(FRperiods(ii,:));
end

% Licks
idx= find(BehaviorAnimals==(cell2mat(AATC_PFC(i).Animal)-7)&BehaviorSession==cell2mat(AATC_PFC(i).TrainingDay))
Session=AATCpfcAllBehavior(idx).Events;
if Session{1}>5
pre=1*1000-1;
post=5*1000;
bin=100;
Fig=0;
[RewardMatrix,unRewardMatrix]=HPCpaper_Pre_BehaviorExample(Session,pre,post,bin,Fig);
[Learned,LearnedP]=ttest2(mean(RewardMatrix(:,pre+2000:pre+2975),2),mean(unRewardMatrix(:,pre+2000:pre+2975),2));
LickRate=mean(sum(RewardMatrix(:,3000:4000),2));
LickRateReward=mean(sum(RewardMatrix(:,4000:5000),2));
else
    Learned=NaN;
    LearnedP=NaN;
    LickRate=NaN;
    LickRateReward=NaN;
end

%% Counters
% Main Psths
AATC_Sua_Psth=cat(2,AATC_Sua_Psth,PsthSua);
%     
TrgDayCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATC_PFC(i).TrainingDay);
AnimalCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATC_PFC(i).Animal);
SessionCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATC_PFC(i).Session);
LearnedCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=Learned;
LearnedPCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LearnedP;
CellPerSesCounter(CellCounter+1:CellCounter+size(PsthSua,2))=1:size(PsthSua,2);
% 
ExpCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATC_PFC(i).Exp);

LickRateCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LickRate;
LickRateRewardCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LickRateReward;

EvokedSigCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=EvokedSig;
TraceSigCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=TraceSig;

AnovaCounter(CellCounter+1:CellCounter+size(PsthSua,2),:)=Anovap;
TaskModCounter(CellCounter+1:CellCounter+size(PsthSua,2))=TaskMod;
FRperiodsCounter=cat(1,FRperiodsCounter,FRperiods);

%%Waveform Extraction
wfsS=cell2mat(AATC_PFC(i).Waveforms.wfs);
ShankS=cell2mat(AATC_PFC(i).Area);
includeWFS=cell2mat(AATC_PFC(i).Waveforms.include);
SpikeV2PS=cell2mat(AATC_PFC(i).Waveforms.SpikeV2P);
SpikeWidthS=cell2mat(AATC_PFC(i).Waveforms.SpikeWidth);
DepthS=cell2mat(AATC_PFC(i).Waveforms.Depth);

%%Waveform Counters
wfs=cat(1,wfs,wfsS(ShankS==1,:));
includeWF=cat(2,includeWF,includeWFS(ShankS==1));
SpikeV2P=cat(1,SpikeV2P,SpikeV2PS(ShankS==1));
spikeWidth=cat(2,spikeWidth,SpikeWidthS(ShankS==1));
Depth=cat(2,Depth,DepthS(ShankS==1)) ; 

CellCounter=CellCounter+size(GoodClusters,1);
clearvars FRperiods PsthSua EventMatrix EventType wfsS includeWFS LickRate SpikeV2PS DepthS Anovap TaskMod ShankS SpikeTS ClusterID GoodClusters RipTS EventType RecordingLength EvokedSig TraceSig

end

clearvars -except  AATC_Sua_Psth TrgDayCounter AnimalCounter LearnedCounter LearnedPCounter...
                    CellPerSesCounter SessionCounter...
                    wfs SpikeV2P spikeWidth Depth includeWF ...
                    AnovaCounter TaskModCounter ...
                    FRperiodsCounter ExpCounter LickRateCounter  LickRateRewardCounter...
                    EvokedSigCounter TraceSigCounter CellCounter CellCounter

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
   Post=4;
   Bin=1;
   
   [PsthSua,EventMatrix,EventType]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,Pre,Post,Bin);
    PFCCells=find(cell2mat(AATC_HPC_PFC(i).Area)==2);
   for ii=1:length(PFCCells)
   EvokedSig(ii)=ranksum(squeeze(mean(EventMatrix(EventType==1,Pre*30000:(Pre+1.35)*30000,ii),2)),...
   squeeze(mean(EventMatrix(EventType==2,Pre*30000:(Pre+1.35)*30000,ii),2))) ;
   TraceSig(ii)=ranksum(squeeze(mean(EventMatrix(EventType==1,(Pre+2)*30000:(Pre+3)*30000,ii),2)),...
   squeeze(mean(EventMatrix(EventType==2,(Pre+2)*30000:(Pre+3)*30000,ii),2))) ;
   end
%% Anova  
        for ii=1:length(PFCCells)
        Baseline=squeeze(mean(EventMatrix(:,1:Pre*30000,ii),2));
        Evoked=squeeze(mean(EventMatrix(:,Pre*30000:Pre*30000+15000,ii),2));
        Sustained=squeeze(mean(EventMatrix(:,Pre*30000+15000:Pre*30000+60000,ii),2));
        Trace=squeeze(mean(EventMatrix(:,Pre*30000+60000:Pre*30000+90000,ii),2));
        Reward=squeeze(mean(EventMatrix(:,Pre*30000+90000:Pre*30000+120000,ii),2));
        Data=cat(1,Baseline,Evoked,Sustained,Trace,Reward);
        TimeGroup=repelem(1:5,size(EventMatrix,1));
        EventTypeGroup=repmat(EventType,5,1);
        [p,tbl,stats] = anovan(Data,{TimeGroup',EventTypeGroup},'model','interaction','display','off');
        Anovap(ii,:)=p;
        if min(Anovap)<0.05
        TaskMod(ii)=1;
        else
        TaskMod(ii)=0;
        end
        end
   
%% Firing Rate
bins=10;
sets=floor(AATC_HPC_PFC(i).RecordingLength/10); 
for ii=1:length(GoodClusters)   
TScell(1:RecordingLength)=0;
TScell(SpikeTS(find(ClusterID==GoodClusters(ii))))=1;
Periods=reshape(TScell(1:sets*(bins)),(bins),sets);
FRperiods(ii,:)=(sum(Periods,2)/(sets/30000))';
FRmean(ii)=mean(FRperiods(ii,:));
FRstd(ii)=std(FRperiods(ii,:));
end
% Licks
idx= find(BehaviorAnimals==(cell2mat(AATC_HPC_PFC(i).Animal)-13)&BehaviorSession==cell2mat(AATC_HPC_PFC(i).TrainingDay))
Session=AATChpcpfcAllBehavior(idx).Events;
pre=1*30000-1;
post=5*30000;
bin=100;
Fig=0;
[RewardMatrix,unRewardMatrix]=HPCpaper_Pre_BehaviorExample(Session,pre,post,bin,Fig);
[Learned,LearnedP]=ttest2(mean(RewardMatrix(:,pre+2000:pre+2975),2),mean(unRewardMatrix(:,pre+2000:pre+2975),2));
LickRate=mean(sum(RewardMatrix(:,3000:4000),2));
LickRateReward=mean(sum(RewardMatrix(:,4000:5000),2));

%% Counters
% % Main Psths
AATC_Sua_Psth=cat(2,AATC_Sua_Psth,PsthSua);
    
TrgDayCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATC_HPC_PFC(i).TrainingDay);
AnimalCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATC_HPC_PFC(i).Animal);
SessionCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATC_HPC_PFC(i).Session);
LearnedCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=Learned;
LearnedPCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LearnedP;
CellPerSesCounter(CellCounter+1:CellCounter+size(PsthSua,2))=1:size(PsthSua,2);

ExpCounter(CellCounter+1:CellCounter+size(PsthSua,2))=cell2mat(AATC_HPC_PFC(i).Exp);

LickRateCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LickRate;
LickRateRewardCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=LickRateReward;
EvokedSigCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=EvokedSig;
TraceSigCounter(CellCounter+1:CellCounter+size(GoodClusters,1))=TraceSig;

AnovaCounter(CellCounter+1:CellCounter+size(PsthSua,2),:)=Anovap;
TaskModCounter(CellCounter+1:CellCounter+size(PsthSua,2))=TaskMod;
FRperiodsCounter=cat(1,FRperiodsCounter,FRperiods);

%%Waveform Extraction
wfsS=cell2mat(AATC_HPC_PFC(i).Waveforms.wfs);
ShankS=cell2mat(AATC_HPC_PFC(i).Area);
includeWFS=cell2mat(AATC_HPC_PFC(i).Waveforms.include);
SpikeV2PS=cell2mat(AATC_HPC_PFC(i).Waveforms.SpikeV2P);
SpikeWidthS=cell2mat(AATC_HPC_PFC(i).Waveforms.SpikeWidth);
DepthS=cell2mat(AATC_HPC_PFC(i).Waveforms.Depth);

%%Waveform Counters
wfs=cat(1,wfs,wfsS(ShankS==1,:));
includeWF=cat(2,includeWF,includeWFS(ShankS==1));
SpikeV2P=cat(1,SpikeV2P,SpikeV2PS(ShankS==1));
spikeWidth=cat(2,spikeWidth,SpikeWidthS(ShankS==1));
Depth=cat(2,Depth,DepthS(ShankS==1)) ; 

CellCounter=CellCounter+size(GoodClusters,1);
clearvars FRperiods PsthSua EventMatrix EventType wfsS includeWFS LickRate SpikeV2PS DepthS Anovap TaskMod ShankS EvokedSig TraceSig

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

AATC_PFC_Sua_Psth=AATC_Sua_Psth;
%% Save
cd('T:\jan\Collabo Data\PFCpaperPreProcessed')
save AATC_PFC_Sua_Psth_1ms AATC_PFC_Sua_Psth TrgDayCounter AnimalCounter ...
                    CellPerSesCounter SessionCounter...
                    wfs SpikeV2P spikeWidth neuronTypeGmm Depth includeWF ...
                    AnovaCounter TaskModCounter ...
                    FRperiodsCounter ExpCounter neuronType ...
                    LearnedCounter LickRateCounter LickRateRewardCounter LearnedPCounter...
                    EvokedSigCounter TraceSigCounter


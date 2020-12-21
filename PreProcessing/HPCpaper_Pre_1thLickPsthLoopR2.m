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
HPCcounter=1;
for i=1:length(AATCR2)
   SpikeTS=AATCR2(i).Spikes{1,1};
   ClusterID=AATCR2(i).Spikes{2,1};
   GoodClusters=unique(ClusterID);
   
  EventTypeAll=AATCR2(i).Events{1,1}{2,1};
   EventTSAll=AATCR2(i).Events{1,1}{1,1};
   
   EventType=EventTypeAll(ismember(EventTypeAll,[1,2,6]));
   EventTS=EventTSAll(ismember(EventTypeAll,[1,2,6]));
   

   
LickTS=EventTSAll(ismember(EventTypeAll,[3]));;
CSTS=EventTSAll(ismember(EventTypeAll,[1,6]));
TimeFromLastLick=5*30000;
TimeFromLastCS=13*30000;
NextSecs=1*30000;



for ii=1:length(LickTS)

    if length(intersect(LickTS(ii)-TimeFromLastLick:LickTS(ii)-1,LickTS))==0&length(intersect(LickTS(ii)-TimeFromLastCS:LickTS(ii)-1,CSTS))==0 ...
       & length(intersect(LickTS(ii):LickTS(ii)+NextSecs,LickTS))>=3;
        LickBoutThres(ii)=1;
    else
        LickBoutThres(ii)=0;
    end
end

% figure()
EventTS=LickTS(LickBoutThres==1);
% EventType(1:length(EventTS))=1;

RecordingLength=AATCR2(i).RecordingLength;

% if length(EventTS)==0
% LickTrigPsthOutsideTrialHPC(HPCcounter:HPCcounter+length(GoodClusters)-1,:)=NaN(length(GoodClusters),4000); 
% else
% [PsthTrigLick,EventMatrix,EventType1]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,2,2,1);
% Bin=1;
% EventMatrixMeans=squeeze(mean(EventMatrix(EventType1==1,:,:)));
% LickTrig=squeeze(sum(reshape(EventMatrixMeans,Bin*30,size(EventMatrixMeans,1)/(Bin*30),length(GoodClusters)))).*(1000/Bin);
% LickTrigPsthOutsideTrialHPC(HPCcounter:HPCcounter+size(PsthTrigLick,2)-1,1,:)=LickTrig';
% 
% EventMatrixMeans=squeeze(mean(EventMatrix(EventType1==2,:,:)));
% LickTrig=squeeze(sum(reshape(EventMatrixMeans,Bin*30,size(EventMatrixMeans,1)/(Bin*30),length(GoodClusters)))).*(1000/Bin);
% LickTrigPsthOutsideTrialHPC(HPCcounter:HPCcounter+size(PsthTrigLick,2)-1,2,:)=LickTrig';
% 
% EventMatrixMeans=squeeze(mean(EventMatrix(EventType1==6,:,:)));
% LickTrig=squeeze(sum(reshape(EventMatrixMeans,Bin*30,size(EventMatrixMeans,1)/(Bin*30),length(GoodClusters)))).*(1000/Bin);
% LickTrigPsthOutsideTrialHPC(HPCcounter:HPCcounter+size(PsthTrigLick,2)-1,3,:)=LickTrig';
% 
% end
% 
clearvars EventType

TimeFromLastLick=5*30000;
TimeFromLastCS=3*30000;
NextSecs=1*30000;
for ii=1:length(LickTS)

    if length(intersect(LickTS(ii)-TimeFromLastLick:LickTS(ii)-1,LickTS))==0&length(intersect(LickTS(ii)-TimeFromLastCS:LickTS(ii)-1,CSTS))==1 ...
       & length(intersect(LickTS(ii):LickTS(ii)+NextSecs,LickTS))>=3;
        LickBoutThres(ii)=1;
    else
        LickBoutThres(ii)=0;
    end
end

% figure()
EventTS=LickTS(LickBoutThres==1);
EventType(1:length(EventTS))=1;


if length(EventTS)==0
LickTrigPsthInsideTrial(HPCcounter:HPCcounter+length(GoodClusters)-1,:)=NaN(length(GoodClusters),4000); 
else
[PsthTrigLick,EventMatrix,EventType1]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,2,2,1);

LickTrigPsthInsideTrial(HPCcounter:HPCcounter+size(PsthTrigLick,2)-1,:)=PsthTrigLick';


end

HPCcounter=HPCcounter+size(GoodClusters,1);
 
clearvars PsthTrigLick EventMatrix EventTypeAll EventTSAll EventType EventTS GoodClusters GoodClusters LickBoutThres EventMatrix EventType1 EventMatrixMeans

end
save HPCPFC_1thLicksInNoutOfTrialsR2 LickTrigPsthInsideTrial


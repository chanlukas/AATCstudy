clear all

%% Load Data from HPConly Experiment
cd('T:\jan\Collabo Data')
load AATCstrucHPConly.mat


%%Exp1
HPCcounter=1;
for i=1:length(AATC_HPC)

   SpikeTS=AATC_HPC(i).Spikes{1,1};
   ClusterID=AATC_HPC(i).Spikes{2,1};
   GoodClusters=unique(ClusterID);
   
   EventTypeAll=AATC_HPC(i).Events{2,1};
   EventTSAll=AATC_HPC(i).Events{1,1};
   

   
LickTS=EventTSAll(ismember(EventTypeAll,[3]));;
CSTS=EventTSAll(ismember(EventTypeAll,[1]));
TimeFromLastLick=5*30000;
TimeFromLastCS=13*30000;
NextSecs=1*30000;

clearvars EventType EventTS 

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
EventType(1:length(EventTS))=1;

RecordingLength=AATC_HPC(i).RecordingLength;

if length(EventTS)==0
LickTrigPsthOutsideTrialHPC(HPCcounter:HPCcounter+length(GoodClusters)-1,:)=NaN(length(GoodClusters),4000); 
else
[PsthTrigLick]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,2,2,1);

LickTrigPsthOutsideTrialHPC(HPCcounter:HPCcounter+size(PsthTrigLick,2)-1,:)=PsthTrigLick';
end

clearvars EventType EventTS LickBoutThres

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
LickTrigPsthInsideTrialHPC(HPCcounter:HPCcounter+length(GoodClusters)-1,:)=NaN(length(GoodClusters),4000); 
else
[PsthTrigLick]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,2,2,1);
LickTrigPsthInsideTrialHPC(HPCcounter:HPCcounter+size(PsthTrigLick,2)-1,:)=PsthTrigLick';
end

HPCcounter=HPCcounter+size(GoodClusters,1);
 
clearvars PsthTrigLick EventMatrix EventTypeAll EventTSAll EventType EventTS GoodClusters GoodClusters LickBoutThres

end

clearvars -except  LickTrigPsthInsideTrialHPC LickTrigPsthOutsideTrialHPC HPCcounter

load AATCstrucPFConly.mat
load AATCpfcAllBehavior.mat

PFCcounter=1;
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
   

   
LickTS=EventTSAll(ismember(EventTypeAll,[3]));;
CSTS=EventTSAll(ismember(EventTypeAll,[1]));
TimeFromLastLick=5*30000;
TimeFromLastCS=13*30000;
NextSecs=1*30000;

clearvars EventType EventTS

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
EventType(1:length(EventTS))=1;
if length(EventTS)==0
LickTrigPsthOutsideTrialPFC(PFCcounter:PFCcounter+length(GoodClusters)-1,:)=NaN(length(GoodClusters),4000); 
else
[PsthTrigLick]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,2,2,1);

LickTrigPsthOutsideTrialPFC(PFCcounter:PFCcounter+size(PsthTrigLick,2)-1,:)=PsthTrigLick';
end

clearvars EventType EventTS LickBoutThres


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
LickTrigPsthInsideTrialPFC(PFCcounter:PFCcounter+length(GoodClusters)-1,:)=NaN(length(GoodClusters),4000); 
else
[PsthTrigLick]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters,EventTS,EventType,RecordingLength,2,2,1);
LickTrigPsthInsideTrialPFC(PFCcounter:PFCcounter+size(PsthTrigLick,2)-1,:)=PsthTrigLick';
end

PFCcounter=PFCcounter+size(GoodClusters,1);

clearvars PsthTrigLick EventMatrix EventTypeAll EventTSAll EventType EventTS GoodClusters GoodClusters LickBoutThres

end

clearvars -except LickTrigPsthInsideTrialHPC LickTrigPsthOutsideTrialHPC LickTrigPsthInsideTrialPFC LickTrigPsthOutsideTrialPFC HPCcounter PFCcounter
load  AATCstrucHPCPFC.mat
load  AATChpcpfcAllBehavior.mat


for i=1:length(AATC_HPC_PFC)
    i
   SpikeTS=AATC_HPC_PFC(i).Spikes{1,1};
   ClusterID=AATC_HPC_PFC(i).Spikes{2,1};
   ShankS=cell2mat(AATC_HPC_PFC(i).Area);
   GoodClusters=unique(ClusterID);
   GoodClusters1=GoodClusters(ShankS==1);
   GoodClusters2=GoodClusters(ShankS==2);

   EventTypeAll=AATC_HPC_PFC(i).Events{2,1};
   EventTSAll=AATC_HPC_PFC(i).Events{1,1};
   
   EventType=EventTypeAll(ismember(EventTypeAll,[1,2]));
   EventTS=EventTSAll(ismember(EventTypeAll,[1,2]));
   
   RecordingLength=AATC_HPC_PFC(i).RecordingLength;
   
LickTS=EventTSAll(ismember(EventTypeAll,[3]));
CSTS=EventTSAll(ismember(EventTypeAll,[1]));
TimeFromLastLick=5*30000;
TimeFromLastCS=13*30000;
NextSecs=1*30000;
clearvars EventType EventTS

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
EventType(1:length(EventTS))=1;


[PsthTrigLick]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters1,EventTS,EventType,RecordingLength,2,2,1);

if size(PsthTrigLick,2)==0
LickTrigPsthOutsideTrialHPC(HPCcounter:HPCcounter+length(GoodClusters1)-1,:)=NaN(length(GoodClusters1),4000); 
else
LickTrigPsthOutsideTrialHPC(HPCcounter:HPCcounter+size(PsthTrigLick,2)-1,:)=PsthTrigLick';
end

[PsthTrigLick]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters2,EventTS,EventType,RecordingLength,2,2,1);

if size(PsthTrigLick,2)==0
LickTrigPsthOutsideTrialPFC(PFCcounter:PFCcounter+length(GoodClusters2)-1,:)=NaN(length(GoodClusters2),4000); 
else
LickTrigPsthOutsideTrialPFC(PFCcounter:PFCcounter+size(PsthTrigLick,2)-1,:)=PsthTrigLick';
end


TimeFromLastLick=5*30000;
TimeFromLastCS=3*30000;
NextSecs=1*30000;

clearvars EventType EventTS LickBoutThres

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
LickTrigPsthInsideTrialHPC(HPCcounter:HPCcounter+length(GoodClusters1)-1,:)=NaN(length(GoodClusters1),4000); 
else
[PsthTrigLick]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters1,EventTS,EventType,RecordingLength,2,2,1);
LickTrigPsthInsideTrialHPC(HPCcounter:HPCcounter+size(PsthTrigLick,2)-1,:)=PsthTrigLick';
end

if length(EventTS)==0
LickTrigPsthInsideTrialPFC(PFCcounter:PFCcounter+length(GoodClusters2)-1,:)=NaN(length(GoodClusters2),4000); 
else
[PsthTrigLick]=HPCpaper_Pre_SuaPsth(SpikeTS,ClusterID,GoodClusters2,EventTS,EventType,RecordingLength,2,2,1);
LickTrigPsthInsideTrialPFC(PFCcounter:PFCcounter+size(PsthTrigLick,2)-1,:)=PsthTrigLick';
end


HPCcounter=HPCcounter+size(GoodClusters1,1);
PFCcounter=PFCcounter+size(GoodClusters2,1);

clearvars FRperiods PsthSua EventMatrix EventType wfsS includeWFS LickRate SpikeV2PS DepthS Anovap TaskMod ShankS EvokedSig TraceSig LickBoutThres

end

%% Save
cd('T:\jan\Collabo Data\PFCpaperPreProcessed')

save HPCPFC_1thLicksInNoutOfTrials LickTrigPsthInsideTrialHPC LickTrigPsthOutsideTrialHPC LickTrigPsthInsideTrialPFC LickTrigPsthOutsideTrialPFC


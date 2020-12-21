
function [ResponseRates]=HPCpaper_BehaviorLoop2R(Events,Animals,Sessions)

%% Prepare Input Variables

animals=cell2mat(Animals);
sessions=cell2mat(Sessions);
%% Loop over animals and sessions per animal
for i= 1:length(unique(animals))
    
for ii=1:length(sessions(animals==(min(animals)-1)+i))
    if i==1&ii==7  %Skip1 session with no behavior data
        continue
    end
%% define the codes for CSstimulus and Lick Response indices, timestamps   
CSstim=[1,2,6];
Response=3;
SessionIdx=find(animals==(min(animals)-1)+i);
EventType=Events{1,SessionIdx(ii)}{2,1};
Time=Events{1,SessionIdx(ii)}{1,1};

%%define binsize and trial length
Bin=100;
pre=1000;
post=5000;

%% sort CS+ and CS- timestamps  

EventIdx=[];
EventTS=[];
for j=1:length(CSstim)
   EventIdx1(1:length(find(EventType==CSstim(j))))=CSstim(j);
   EventIdx=horzcat(EventIdx,EventIdx1);
   EventTS1=Time(find(EventType==CSstim(j)));
   EventTS=vertcat(EventTS1,EventTS);   
   clearvars EventIdx1 EventTS1
end

%% Make Licktimestamps into continous Lick vector
ResponseTS(1:Time(end))=0;
idxq=find(EventType==Response);
ResponseTS(Time(idxq))=1;

%% compute peri stimulus lick matrix

EventResponseMatrix=[];  
EventIdxTimeOK=[];
for j=1:length(EventTS)
    if (EventTS(j)-(pre-1))>0&&(EventTS(j)+post<length(ResponseTS))      %filtering out events for which the pre post cut off would violate recording bounds    
        EventResponseMatrix=cat(1,EventResponseMatrix,ResponseTS(EventTS(j)-(pre-1):EventTS(j)+post));        
        EventIdxTimeOK(j)=1;
    else
        EventIdxTimeOK(j)=0;
    end
end


%% Compute Significance
% EventIdx1=EventIdx(EventIdxTimeOK==1);
% [Learned{i,ii},LearnedP]=ttest2(nanmean(EventResponseMatrix(EventIdx1==CSstim(1),pre+2000:pre+2975),2),nanmean(EventResponseMatrix(find(EventIdx1==CSstim(2)),pre+2000:pre+2975),2));   

%% Average and Bin Lick Responses and convert into HZ
for j=1:size(CSstim,2)
meanEventResponses=mean(EventResponseMatrix(EventIdx(EventIdxTimeOK==1)==CSstim(j),:)',2);
binMeanEventResponsesHZ(j,:)=sum(reshape(meanEventResponses,Bin,size(meanEventResponses,1)/Bin))*(1000/Bin);
clearvars meanEventResponses 
end

ResponseRatesSingleSession(ii,:,:)=binMeanEventResponsesHZ;

clearvars EventIdxTimeOK
end
ResponseRates{i}=ResponseRatesSingleSession;
clearvars ResponseRatesSingleSession
end
end



function [RewardMatrix,RewardMatrix2,unRewardMatrix]=HPCpaper_BehaviorExampleR2(Session,pre,post,bin)

EventType=Session{2,1};
Time=Session{1,1};

RewardTS1=Time(find(EventType==1));
unRewardTS=Time(find(EventType==2));
RewardTS2=Time(find(EventType==6));

LickTS(1:Time(end))=0;
LickTSraw=(Time(find(EventType==3)));
LickTSraw(isnan(LickTSraw))=1;  
LickTS(LickTSraw)=1;


RewardMatrix=[];  
for i=1:length(RewardTS1)
    if (RewardTS1(i)-pre)>0&&(RewardTS1(i)+post<length(LickTS))      %filtering out events for which the pre post cut off would violate recording bounds    
        RewardMatrix=cat(1,RewardMatrix,LickTS(RewardTS1(i)-pre:RewardTS1(i)+post));        
    end    
end
RewardMatrix(:,4995:5095)=RewardMatrix(:,4895:4995);

RewardMatrix2=[];  
for i=1:length(RewardTS2)
    if (RewardTS2(i)-pre)>0&&(RewardTS2(i)+post<length(LickTS))      %filtering out events for which the pre post cut off would violate recording bounds    
        RewardMatrix2=cat(1,RewardMatrix2,LickTS(RewardTS2(i)-pre:RewardTS2(i)+post));        
    end    
end
RewardMatrix2(:,4995:5095)=RewardMatrix2(:,4895:4995);


unRewardMatrix=[];  
for i=1:length(unRewardTS)
    if (unRewardTS(i)-pre)>0&&(unRewardTS(i)+post<length(LickTS))      %filtering out events for which the pre post cut off would violate recording bounds    
        unRewardMatrix=cat(1,unRewardMatrix,LickTS(unRewardTS(i)-pre:unRewardTS(i)+post));        
    end    
end



meanRewardLicks=mean(RewardMatrix(:,:)',2);
binMeanRewardLicks=sum(reshape(meanRewardLicks,bin,size(RewardMatrix,2)/bin))*10;

meanRewardLicks2=mean(RewardMatrix2(:,:)',2);
binMeanRewardLicks2=sum(reshape(meanRewardLicks2,bin,size(RewardMatrix2,2)/bin))*10;

meanUnRewardLicks=mean(unRewardMatrix(:,:)',2);
binMeanUnRewardLicks=sum(reshape(meanUnRewardLicks,bin,size(RewardMatrix,2)/bin))*10;



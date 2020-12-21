

function [RewardMatrix,unRewardMatrix]=HPCpaper_BehaviorExample(Session,pre,post,bin,Fig)

EventType=Session{2,1};
Time=Session{1,1};

RewardTS=Time(find(EventType==1));
unRewardTS=Time(find(EventType==2));

LickTS(1:Time(end))=0;
LickTSraw=(Time(find(EventType==3)));
LickTSraw(isnan(LickTSraw))=1;  
LickTS(LickTSraw)=1;


RewardMatrix=[];  
for i=1:length(RewardTS)
    if (RewardTS(i)-pre)>0&&(RewardTS(i)+post<length(LickTS))      %filtering out events for which the pre post cut off would violate recording bounds    
        RewardMatrix=cat(1,RewardMatrix,LickTS(RewardTS(i)-pre:RewardTS(i)+post));        
    end    
end


% RewardMatrix(:,4995:5095)=RewardMatrix(:,4895:4995);

unRewardMatrix=[];  
for i=1:length(unRewardTS)
    if (unRewardTS(i)-pre)>0&&(unRewardTS(i)+post<length(LickTS))      %filtering out events for which the pre post cut off would violate recording bounds    
        unRewardMatrix=cat(1,unRewardMatrix,LickTS(unRewardTS(i)-pre:unRewardTS(i)+post));        
    end    
end



meanRewardLicks=mean(RewardMatrix(:,:)',2);
binMeanRewardLicks=sum(reshape(meanRewardLicks,bin,size(RewardMatrix,2)/bin))*10;
meanUnRewardLicks=mean(unRewardMatrix(:,:)',2);
binMeanUnRewardLicks=sum(reshape(meanUnRewardLicks,bin,size(RewardMatrix,2)/bin))*10;

if Fig==1
figure()
% title('Example Session')
subplot(2,2,[3,4])
f=fill([0,0,2000,2000],[0,max(binMeanRewardLicks),max(binMeanRewardLicks),0],[.7,.7,.9])
set(f,'LineStyle','none')
hold on
f1=fill([3000,3000,3015,3015],[0,max(binMeanRewardLicks),max(binMeanRewardLicks),0],[.7,.9,.7])
set(f1,'LineStyle','none')
x=-pre:bin:post;
h1=plot(x,binMeanRewardLicks,'b')
set(h1,'LineWidth',3)
hold on
h2=plot(x,binMeanUnRewardLicks,'r')
set(h2,'LineWidth',3)
set(gca,'FontSize',14)
xlabel('Time (s)')
ylabel('Licks (hz)')
axis tight
box off
set(gca,'FontSize',20)

TS1=RewardTS;
TS2=find(LickTS==1);
[Trial,Timestamps]=Rasterplot(TS1,TS2,pre,post)
find(Timestamps>3000&Timestamps<3010)
Timestamps(find(Timestamps>3000&Timestamps<3050))=NaN;
subplot(2,2,1)

f=fill([0,0,2000,2000],[0,size(Trial,2),size(Trial,2),0],[.7,.7,.9])
set(f,'LineStyle','none')
hold on
f1=fill([3000,3000,3050,3050],[0,size(Trial,2),size(Trial,2),0],[.7,.9,.7])
set(f1,'LineStyle','none')
hold on
h=plot(Timestamps,Trial,'s','markersize',2,'markerfacecolor',[0 0 1],'markeredgecolor',[0 0 1]);
axis([-pre post 0 size(TS1,1)])
xlabel('Time (s)')
ylabel('Trials')
set(gca,'FontSize',20)
box off

TS1=unRewardTS;
TS2=find(LickTS==1);
[Trial1,Timestamps1]=Rasterplot(TS1,TS2,pre,post)

subplot(2,2,2)
f=fill([0,0,2000,2000],[0,size(TS1,1),size(TS1,1),0],[.7,.7,.9])
set(f,'LineStyle','none')
hold on
f1=fill([3000,3000,3050,3050],[0,size(Trial,2),size(Trial,2),0],[.7,.9,.7])
set(f1,'LineStyle','none')
hold on
h=plot(Timestamps1,Trial1,'s','markersize',2,'markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0]);
axis([-pre post 0 size(TS1,1)])
xlabel('Time (s)]')
ylabel('Trials')
set(gca,'FontSize',20)
box off
end



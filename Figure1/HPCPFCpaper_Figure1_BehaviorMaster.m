clear all

%% Load and Combine All Event Data from HPConly Experiment
cd('T:\jan\Collabo Data')
load AATChpcAllBehavior.mat
Events={AATChpcAllBehavior.Events};
Animals={AATChpcAllBehavior.Animal};
Sessions={AATChpcAllBehavior.Session};

[ResponseRates]=HPCpaper_Pre_BehaviorLoop(Events,Animals,Sessions);

clearvars -except ResponseRates

load AATChpcpfcAllBehavior.mat
Events={AATChpcpfcAllBehavior.Events};
Animals={AATChpcpfcAllBehavior.Animal};
Sessions={AATChpcpfcAllBehavior.Session};
[ResponseRates1]=HPCpaper_Pre_BehaviorLoop(Events,Animals,Sessions);

clearvars -except ResponseRates ResponseRates1

load  AATCpfcAllBehavior.mat
Events={AATCpfcAllBehavior.Events};
Animals={AATCpfcAllBehavior.Animal};
Sessions={AATCpfcAllBehavior.Session};
[ResponseRates2]=HPCpaper_Pre_BehaviorLoop(Events,Animals,Sessions);


clearvars -except ResponseRates ResponseRates1 ResponseRates2
%% Plot Average Lick Responses

ResponseRates=[ResponseRates,ResponseRates1,ResponseRates2];
figure()
for i=1:size(ResponseRates,2)

ResponseRates1=squeeze(ResponseRates{i});
ResponseRatesTraces(i,:,:)=squeeze(nanmean(ResponseRates1(1:7,:,30:40),3))-squeeze(nanmean(ResponseRates1(1:7,:,1:9),3));


h=plot(squeeze(ResponseRatesTraces(i,:,1)),'b')
set(h,'LineWidth',3)
set(gca,'FontSize',20)
hold on
h=plot(squeeze(ResponseRatesTraces(i,:,2)),'r')
set(h,'LineWidth',3)
set(gca,'FontSize',20)
xlabel('Session')
ylabel('LickRate Change Trace Interval [Hz]')
title('Session Overview')
axis tight
box off
hold on
end

figure()
pp(1,:)=nanmean(ResponseRatesTraces(:,:,2));
errBar=repmat(std(ResponseRatesTraces(:,:,2))/sqrt(size(ResponseRatesTraces,1)),2,1);
shadePlot2(1:7,pp(1,:),errBar,[0 0 1],.3)
hold on
box off
axis tight

pp(2,:)=nanmean(ResponseRatesTraces(:,:,1));
errBar=std(ResponseRatesTraces(:,:,1))/sqrt(size(ResponseRatesTraces,1));
shadePlot2(1:7,pp(2,:),errBar',[1 0 0],.3)
hold on
ylabel('{\Delta} LickRate Trace Period [hz]')
xlabel('Training Days')

for i=1:7
  [p(i),h]=ranksum(ResponseRatesTraces(:,i,2),ResponseRatesTraces(:,i,1))
end

box off
axis tight

clear all
load AATChpcpfcAllBehavior.mat
Events={AATChpcpfcAllBehavior.Events};
Animals={AATChpcpfcAllBehavior.Animal};
Sessions={AATChpcpfcAllBehavior.Session};

PreSession=Events{1,find(cell2mat(Animals)==3&cell2mat(Sessions)==1)};
pre=2000;
post=4999;
bin=100;
Fig=1
HPCpaper_Pre_BehaviorExample(PreSession,pre,post,bin,Fig)

PostSession=Events{1,find(cell2mat(Animals)==3&cell2mat(Sessions)==8)};
pre=2000;
post=4999;
bin=100;
Fig=1
HPCpaper_Pre_BehaviorExample(PostSession,pre,post,bin,Fig)




% % Code for reproducing Figure 5
% % Requires Cell-Assembly-Detection package on path (see https://github.com/tortlab/Cell-Assembly-Detection)
% % Uses sub.m auxiliary function to plot.

%%
%% Part 1

clear
load('preprocessed_data', 'data', 'hasPF', 'hasCA1', 'learned', 'ses', 'exp', 'trainingday', 'animal',...
    'varnames', 'srate', 'datafiles')

%%

CA1=1;
PFC=2;

PRE = 1;
STIM = 2;
TRACE = 3;
RWD = 4;
ALL = 5;

timebins = -1:0.001:4;
periods={};
periods{PRE} = find(timebins<0);
periods{STIM} = find(timebins>=0 & timebins<1);
periods{TRACE} = find(timebins>=2 & timebins<3);
periods{RWD} = find(timebins>=3 & timebins<4);
periods{ALL} = find(timebins>=-Inf);


doplot = false;
valid_ses=[hasCA1', hasPF', hasCA1' & hasPF'];
valid_neurons={CA1, PFC, []};
nsub = 20;
binwitdh = 20;
nneurons_all = [];
ntrials_all = [];

area_labels={'CA1', 'PFC', 'CA1 & PFC'};

try
    load(['Assembly_patterns'])
    computeassemblies=0;
catch
    computeassemblies=1;
end

try 
    
    load(['Figure5_data1'], 'area_labels',...
        'CA1','PFC','binwitdh','valid_neurons','valid_ses',...
        'mean_assembly_all','seed')
catch

    mean_assembly_all = {};
    for ises = 1:size(valid_ses,1) % find(hasPF & ~hasCA1)% [1:60 92:size(valid_ses,1)]
        seed = 654299090;
        clear mean_assembly
    %     return
    %     rng('default') 
        rng(seed)
        
        animal(ises) = data(ises).Animal{1};
        spktimes = double(data(ises).Spikes{1})/srate;
        spkid = data(ises).Spikes{2};
        area = data(ises).Area{1};
        etime = data(ises).Events{1}/srate;
        etype = data(ises).Events{2};
        spkid = spkid-min(spkid)+1;
        classes = unique(spkid);

        if isempty(data(ises).RipTS)
            continue
        end
        riptime = double(data(ises).RipTS{1})/srate;

        area_aux =[];
        area_aux(classes) =area;


        nneurons=[];
        for iarea=1:length(area_labels)-1
            nneurons(iarea) = sum(area_aux(unique(spkid))==valid_neurons{iarea});
        end
        if sum(nneurons==0)
            nneurons(iarea+1) = 0;
            
        else
            nneurons(iarea+1)= sum(nneurons(1:2));
        end


        idx1 = find(etype == 1);
        trialtype = 1*ones(length(idx1),1);
        idx2 = find(etype == 2);
        trialtype = [trialtype; 2*ones(length(idx2),1)];
        idx = [idx1; idx2];
        trialtime = etime(idx);
        trialtime(:,2) = etime(idx)+3;

        timebins = -1:0.001:4;


        %%
        for iarea = find(nneurons~=0)

            disp(area_labels{iarea})

            nsubsample=1;

            for isample = 1:nsubsample
                area_=area;
                area_both=area;

                valid = 0;
                timebins = -1:0.001:4;

                for isource=1:length(periods)
                %%
                timebins = -1:0.001:4;
                psth=zeros(length(classes),size(trialtime,1),length(timebins));
                psth_smooth=psth;
                for itrial = 1:size(trialtime,1)
                    valid = spktimes>=trialtime(itrial,1)-1 & spktimes<trialtime(itrial,2)+1;

                    if ~isempty(valid_neurons{iarea})
                        area_aux =[];
                        area_aux(classes) =area_;
                        valid = valid & (area_aux(spkid)==valid_neurons{iarea})';
                    else
                        area_aux =[];
                        area_aux(classes) =area_both;
                        valid = valid & (area_aux(spkid)>0)';
                    end
                    spktimes_ = spktimes(valid);
                    spktimes_ = spktimes_-trialtime(itrial,1);

                    spkid_ = spkid(valid);

                    for ineuron=1:length(classes)
                        idx = find(spkid_==classes(ineuron));
                        if isempty(idx)
                            continue
                        end
                        count = histc(spktimes_(idx),timebins);
                        psth(ineuron,itrial,:)=count;
                        psth_smooth(ineuron,itrial,:)=conv(count,rectwin(binwitdh),'same');
                    end
                end
                %%
                
                if computeassemblies==1
                    M = psth_smooth(:,:,periods{isource});
                    M = M(:,:,1:20:end);
                    M = reshape(M,size(M,1),[]);
                    P = assembly_patterns(M);
                else
                    P=patterns_all{iarea,ises,isource};
                end


                psth_proj=zeros(size(P,2),size(psth_smooth,2),length(timebins));
                parfor itrial=1:size(psth_smooth,2)
                    psth_proj(:,itrial,:) = assembly_activity(P,squeeze(psth_smooth(:,itrial,:)));
                end            

                if isempty(P)
                    continue
                end

                mean_assembly_stim = squeeze(nanmean(psth_proj,2));
                sem_assembly_stim = squeeze(nanstd(psth_proj,0,2)./...
                    sqrt(sum(~isnan(psth_proj),2)));
                mean_assembly_stim1 = squeeze(nanmean(psth_proj(:,trialtype==1,:),2));
                sem_assembly_stim1 = squeeze(nanstd(psth_proj(:,trialtype==1,:),0,2)./...
                    sqrt(sum(~isnan(psth_proj(:,trialtype==1,:)),2)));
                mean_assembly_stim2 = squeeze(nanmean(psth_proj(:,trialtype==2,:),2));
                sem_assembly_stim2 = squeeze(nanstd(psth_proj(:,trialtype==2,:),0,2)./...
                    sqrt(sum(~isnan(psth_proj(:,trialtype==2,:)),2)));

                %%
                timebins = -1:0.001:1;
                licktime = etime(etype==3);
                licktime(:,2) = licktime+1;
                clear psth psth_smooth
                psth = zeros(length(classes),size(licktime,1),length(timebins));
                psth_smooth = zeros(length(classes),size(licktime,1),length(timebins));

                for ilick = 1:size(licktime,1)
                    valid = spktimes>=licktime(ilick,1)-1 & spktimes<licktime(ilick,2);

                    if ~isempty(valid_neurons{iarea})
                        area_aux =[];
                        area_aux(classes) =area_;
                        valid = valid & (area_aux(spkid)==valid_neurons{iarea})';
                    else
                        area_aux =[];
                        area_aux(classes) =area_both;
                        valid = valid & (area_aux(spkid)>0)';
                    end
                    spktimes_ = spktimes(valid);
                    spktimes_ = spktimes_-licktime(ilick,1);

                    spkid_ = spkid(valid);

                    for ineuron=1:length(classes)
                        idx = find(spkid_==classes(ineuron));
                        if isempty(idx)
                            continue
                        end
                        count = histc(spktimes_(idx),timebins);
                        psth(ineuron,ilick,:)=count;
                        psth_smooth(ineuron,ilick,:)=conv(count,rectwin(binwitdh),'same');
                    end
                    %%
                end
    %             toc
                %%
                psth_proj_lick=zeros(size(P,2),size(psth_smooth,2),length(timebins));
                parfor ilick=1:size(psth_smooth,2)
                    psth_proj_lick(:,ilick,:) = assembly_activity(P,squeeze(psth_smooth(:,ilick,:)));
                end            


                mean_assembly_lick = squeeze(nanmean(psth_proj_lick,2));
                sem_assembly_lick = squeeze(nanstd(psth_proj_lick,0,2)./sqrt(sum(~isnan(psth_proj_lick),2)));
                %%

                timebins = -1:0.001:1;
                riptime = double(data(ises).RipTS{1})/srate;
                riptime=riptime(:);
                riptime(:,2) = riptime+1;
                clear psth psth_smooth
                psth = zeros(length(classes),size(riptime,1),length(timebins));
                psth_smooth = zeros(length(classes),size(riptime,1),length(timebins));
                for irip = 1:size(riptime,1)
                    valid = spktimes>=riptime(irip,1)-1 & spktimes<riptime(irip,2);

                    if ~isempty(valid_neurons{iarea})
                        area_aux =[];
                        area_aux(classes) =area_;
                        valid = valid & (area_aux(spkid)==valid_neurons{iarea})';
                    else
                        area_aux =[];
                        area_aux(classes) =area_both;
                        valid = valid & (area_aux(spkid)>0)';
                    end
                    spktimes_ = spktimes(valid);
                    spktimes_ = spktimes_-riptime(irip,1);

                    spkid_ = spkid(valid);

                    for ineuron=1:length(classes)
                        idx = find(spkid_==classes(ineuron));
                        if isempty(idx)
                            continue
                        end
                        count = histc(spktimes_(idx),timebins);
                        psth(ineuron,irip,:)=count;
                        psth_smooth(ineuron,irip,:)=conv(count,rectwin(binwitdh),'same');
                  
                    end
                    %%
                end
                %%
                psth_proj_rip=zeros(size(P,2),size(psth_smooth,2),length(timebins));
                parfor irip=1:size(psth_smooth,2)
                    psth_proj_rip(:,irip,:) = assembly_activity(P,squeeze(psth_smooth(:,irip,:)));
                end


                mean_assembly_rip = squeeze(nanmean(psth_proj_rip,2));
                sem_assembly_rip = squeeze(nanstd(psth_proj_rip,0,2)./sqrt(sum(~isnan(psth_proj_rip),2)));

                %%

                mean_assembly{1,1} = mean_assembly_stim1;
                mean_assembly{2,1} = mean_assembly_stim2;
                mean_assembly{3,1} = mean_assembly_stim;
                mean_assembly{4,1} = mean_assembly_lick;
                mean_assembly{5,1} = mean_assembly_rip;

                mean_assembly{1,2} = sem_assembly_stim1;
                mean_assembly{2,2} = sem_assembly_stim2;
                mean_assembly{3,2} = sem_assembly_stim;
                mean_assembly{4,2} = sem_assembly_lick;
                mean_assembly{5,2} = sem_assembly_rip;

                mean_assembly_all{iarea,ises,isource} = mean_assembly;

                end

            if isempty(lessneuron)
                nneurons_all(iarea,ises) = size(psth,1);
            else
                nneurons_all(iarea,ises) = lessneuron;
            end

            nass_all(iarea,ises) = size(P,2);

            ntrial = min(sum(trialtype==1),sum(trialtype==2));
            ntrials_all(iarea,ises) = ntrial;


            end     
        end

        %%
        save(['Figure5_data1'], 'area_labels',...
            'CA1','PFC','binwitdh','valid_neurons','valid_ses',...
            'mean_assembly_all','seed','-v7.3')

    end
end
%%
   
PRE = 1;
STIM = 2;
TRACE = 3;
RWD = 4;
ALL = 5;
source_label = {'PRE','STIM','TRACE','RWD','ALL'};

CSPLUS=1;
CSMINUS=2;
STIM=3;
LICK=4;
RIP = 5;
learn_labels = {'learning' 'overtrained'}

col = [.2 .2 .2;
    .2 .7 .2];

figure(5);clf


sub(7,4,7,1)
hold on
plot(0,nan,'-','linewidth',2,'color',col(1,:))
plot(0,nan,'-','linewidth',2,'color',col(2,:))
% %
aux_stat={};
W_all={};
for iarea = 1:2%3
    
    for isource = ALL
    for ilearn = 1:2
        
        %%
        idx = find(learned(:,1)==ilearn-1 & valid_ses(:,iarea));
        ises = idx;

        m_all_rip=[];
        m_all_stim1=[];
        m_all_stim2=[];
        m_all_lick=[];
        for ises_ = ises'

            mean_assembly = mean_assembly_all{iarea,ises_,isource};

            if isempty(mean_assembly)
                continue
            end
            
            P = patterns_all{iarea,ises_,isource};

            ass_norm=[];
            for iass = 1:size(P,2)
                w=P(:,iass);
                T = w*w';
                T = T-diag(diag(T)); 
                ass_norm(iass) = ((w)'*T*(w));
            end
             
            m = mean_assembly{RIP,1};
            if size(m,2)==1
                m=m';
            end
            
            m = bsxfun(@rdivide, m, ass_norm');
            
            m_all_rip = cat(1,m,m_all_rip);
            if size(m_all_rip,1)>=69
%                 return
            end
            m = mean_assembly{LICK,1};
            if size(m,2)==1
                m=m';
            end
            
            m = bsxfun(@rdivide, m, ass_norm');
            m_all_lick= cat(1,m,m_all_lick);
            
            m = mean_assembly{CSPLUS,1};
            if size(m,2)==1
                m=m';
            end
            
            m = bsxfun(@rdivide, m, ass_norm');
            m_all_stim1= cat(1,m,m_all_stim1);
            
            m = mean_assembly{CSMINUS,1};
            if size(m,2)==1
                m=m';
            end
            
            m = bsxfun(@rdivide, m, ass_norm');
            m_all_stim2= cat(1,m,m_all_stim2);

        end
        
        %%
        
        timebins = -1:0.001:1;
        
        idx_stim = find(timebins>=-2 & timebins<=0);
        idx = find(timebins>=-0.1 & timebins<=0.1);
        
        
        norm = cat(2,m_all_stim1,m_all_stim2);
        m = mean(norm(:,[idx_stim idx_stim+size(m_all_stim1,2)]),2);
        s = std(norm(:,[idx_stim idx_stim+size(m_all_stim1,2)]),0,2);
        
        m_all_stim1 = bsxfun(@rdivide,bsxfun(@minus,m_all_stim1,m),s);
        m_all_stim2 = bsxfun(@rdivide,bsxfun(@minus,m_all_stim2,m),s);
        
        idx = find(timebins>=-0.025 & timebins<=0.025);
        norm = m_all_rip(:,timebins<=-0.05 | timebins>=0.05);
        z_m_all_rip  = bsxfun(@rdivide, bsxfun(@minus, m_all_rip, mean(norm,2)), std(norm,0,2));
        z_m_all_rip = z_m_all_rip(find(~isnan(sum(z_m_all_rip,2))),:);
        
        idx = find(timebins>=-0.05 & timebins<=0.05);
        
        peak = mean(z_m_all_rip(:,idx),2);
        [v, idx] = sort(peak);
        
        sub(7,4,(ilearn-1)*3+(1:3),iarea)
        imagesc(timebins,1:size(m_all_rip,1),z_m_all_rip(idx,:))
        colorbar
        xlim([-0.15 .15])

        if iarea==1
            caxis([-30 200])
        else
            caxis([-10 80])
        end
        
        if ilearn==2
            xlabel('Time from ripple (s)')
        end
        ylabel('Assembly #')
        title([area_labels{iarea} '-' learn_labels{ilearn} ': ' source_label{isource} ])
        
        
        sub(7,4,7,iarea)
        hold on
        
        m=mean((z_m_all_rip(idx,:)')',1);
        s=std((z_m_all_rip(idx,:)')')./sqrt(sum(~isnan((z_m_all_rip(idx,:)')')));
        
        plot(timebins,m,'linewidth',2,'color',col(ilearn,:))
        plot(timebins,m+s,'linewidth',1,'color',col(ilearn,:))
        plot(timebins,m-s,'linewidth',1,'color',col(ilearn,:))
        
        xlim([-0.15 .15])
        xlabel('Time from ripple (s)')
        ylabel({'Mean assembly'; 'reactivation'})
        
        %%
        
        
        
        idx = find(timebins>=-0.025 & timebins<=0.025);
        W = mean(z_m_all_rip(:,idx),2);
        W(isnan(W))=[];
        W_all{iarea,ilearn} = W;
        
        idx = find(timebins>=-0.05 & timebins<=0.05);
        aux_stat{iarea,ilearn} = z_m_all_rip;
        
        
        if iarea==1
            edges = -20:0.2:100;
        else
            edges = -5:0.2:15;
        end
        
        counts = histc(W,edges);
        counts = counts/sum(counts);
        
        sub(5,6,iarea,4:5)
        hold on
        bar(edges,counts,'facecolor','none','edgecolor',col(ilearn,:),'linewidth',2)
        ylabel('PDF')
        title([area_labels{iarea} '-' learn_labels{ilearn} ': ' source_label{isource} ])
        axis tight
        xlabel('Mean Rip. Reactivation')
        
        if iarea==1
            edges = -200:0.2:200;
        else
            edges = -100:0.2:100;
        end
        
        counts = histc(W,edges);
        counts = counts/sum(counts);
        
        sub(5,6,iarea,6)
        hold on
        plot(edges,cumsum(counts),'-','color',col(ilearn,:),'linewidth',2)
        ylabel('CDF')
        title([area_labels{iarea} '-' learn_labels{ilearn} ': ' source_label{isource} ])
        axis tight
        xlabel({'Mean ripple' ; 'reactivation'})
        legend({'Pre' 'Post'},'Location','SouthEast')
    end


    sub(5,6,iarea,6)
    [h, p] = kstest2(W_all{iarea,1},W_all{iarea,2},'tail','larger')
    if iarea==1
        xlim([-25,100])
    else
        xlim([-5,25])
    end
    title([area_labels{iarea} ' - p=' num2str(p)])

    end
    
    if iarea==1
        legend({'Pre', 'Post'},'box','off')
    end
    
    colorbar
    xlim([-0.15 .15])
end


%% Part 2

clear
load('preprocessed_data', 'data', 'hasPF', 'hasCA1', 'learned', 'ses', 'exp', 'trainingday', 'animal',...
    'varnames', 'srate', 'datafiles')

%%

day=[];
animal=[];
for i=1:length(data)
    day(i)= data(i).TrainingDay{1}(1);
    animal(i)= data(i).Animal{1}; 
end

%%
nlick_all={};
for ises = 1:length(data)
    disp(['Performing session: ' num2str(ises)])
    animal(ises) = data(ises).Animal{1};
    spktimes = double(data(ises).Spikes{1})/srate;
    spkid = data(ises).Spikes{2};
    area = data(ises).Area{1};
    etime = data(ises).Events{1}/srate;
    etype = data(ises).Events{2};
    spkid = spkid-min(spkid)+1;
    classes = unique(spkid);

    %%
    idx1 = find(etype == 1);
    trialtype = 1*ones(length(idx1),1);
    idx2 = find(etype == 2);
    trialtype = [trialtype; 2*ones(length(idx2),1)];
    idx = [idx1; idx2];
    trialtime = etime(idx);
    trialtime(:,2) = etime(idx)+3;
    
    valid = 0;

    step=200/1000;
    timebins = -1:step:4;
    binwidth=step;
    preidx = find(timebins>-1 & timebins+binwidth<=0);
    traceidx = find(timebins>2 & timebins+binwidth<=3);
    
    nlick=[];
    nlick_isi=[];
    nlick_first=[];
    for itrial = 1:size(trialtime,1)
        valid = etime>=trialtime(itrial,1)-1 & etime<trialtime(itrial,2)+1;
        valid = valid & etype==3;
%         return
        etime_ = etime(valid);
        etime_ = etime_-trialtime(itrial,1);
        etype_ = etype(valid);
        
        for ibin=1:length(timebins)-1
            idx = find(etime_>timebins(ibin) & etime_<=timebins(ibin+1));
            nlick(itrial,ibin) = length(idx);
            
        end
        idx = find(etime_>0 & etime_<=2);
        nlick_isi(itrial,1) = mean(diff(etime_(idx)));
        
        nlick_first(itrial,1)=NaN;
        if ~isempty(idx)
            nlick_first(itrial,1) = min(etime_(idx));
        end
        
        idx = find(etime_>2 & etime_<=3);
        nlick_isi(itrial,2) = mean(diff(etime_(idx)));
        
        nlick_first(itrial,2)=NaN;
        if ~isempty(idx)
            nlick_first(itrial,2) = min(etime_(idx));
        end
        
        idx = find(etime_>3 & etime_<=4);
        nlick_isi(itrial,3) = mean(diff(etime_(idx)));
        nlick_first(itrial,3)=NaN;
        if ~isempty(idx)
            nlick_first(itrial,3) = min(etime_(idx));
        end
        
        
        %%
    end
    trialtype_all{ises} = trialtype;
    trialtime_all{ises} = trialtime;
    nlick_all{ises} = nlick;
    nlick_isi_all{ises} =nlick_isi;
    nlick_first_all{ises} =nlick_first;
end


idx_pre = find(timebins<-0.1);
idx_stim = find(timebins>=-0.1 & timebins<2);
idx_trace = find(timebins>=2 & timebins<3);
idx_rwd = find(timebins>=3);
%%

CA1=1;
PFC=2;

PRE = 1;
STIM = 2;
TRACE = 3;
RWD = 4;
ALL = 5;


timebins = -1:0.001:4;
periods={};
periods{PRE} = find(timebins<0);
periods{STIM} = find(timebins>=0 & timebins<0.5);
periods{TRACE} = find(timebins>=2 & timebins<3);
periods{RWD} = find(timebins>=3 & timebins<3.5);
periods{ALL} = find(timebins>=-Inf);


doplot = false;
valid_ses=[hasCA1', hasPF', hasCA1' & hasPF'];
valid_neurons={CA1, PFC, []};
nsub = 20;
binwitdh = 20;

area_labels={'CA1', 'PFC', 'CA1 & PFC'};

%%
clear mean_assembly_all 

try
    load(['Figure5_data2'])
catch

    mean_assembly_all={};

    % 0 = use precomputed assemblies. 1 = recompute assemblies
    try
        load(['Assembly_patterns'])
        computeassemblies=0;
    catch
        computeassemblies=1;
    end

    %%
    for ises = 1:size(valid_ses,1) 
        %%% Seed used to compute this analyses
        %seed = 654299090;
        %rng(seed)

        clear mean_assembly

        animal(ises) = data(ises).Animal{1};
        spktimes = double(data(ises).Spikes{1})/srate;
        spkid = data(ises).Spikes{2};
        area = data(ises).Area{1};
        etime = data(ises).Events{1}/srate;
        etype = data(ises).Events{2};
        spkid = spkid-min(spkid)+1;
        classes = unique(spkid);

        if isempty(data(ises).RipTS)
            continue
        end
        riptime = double(data(ises).RipTS{1})/srate;


        area_aux =[];
        area_aux(classes) =area;

        nneurons=[];
        for iarea=1:length(area_labels)-1
            nneurons(iarea) = sum(area_aux(unique(spkid))==valid_neurons{iarea});
        end

        if sum(nneurons==0)
            nneurons(iarea+1) = 0;
        else
            nneurons(iarea+1)= sum(nneurons(1:2));
        end    

        % Finding CS+ and CS- stimulus events
        idx1 = find(etype == 1);
        trialtype = 1*ones(length(idx1),1);
        idx2 = find(etype == 2);
        trialtype = [trialtype; 2*ones(length(idx2),1)];
        idx = [idx1; idx2];
        trialtime = etime(idx);
        trialtime(:,2) = etime(idx)+3;

        [v, trialorder] = sort(trialtime(:,1));

        timebins = -1:0.001:4;


        %%
        for iarea = find(nneurons~=0)

            % No need to subsample in this analysis
            nsubsample=1;
            for isample = 1:nsubsample
                area_=area;
                area_both=area;

                valid = 0;


                for isource=ALL % uses Assembly Patterns computed using all trial period
                timebins = -1:0.001:4;
                psth=zeros(length(classes),size(trialtime,1),length(timebins));
                psth_smooth=psth;
                for itrial = 1:size(trialtime,1)
                    valid = spktimes>=trialtime(itrial,1)-1 & spktimes<trialtime(itrial,2)+1;

                    if ~isempty(valid_neurons{iarea})
                        area_aux =[];
                        area_aux(classes) =area_;
                        valid = valid & (area_aux(spkid)==valid_neurons{iarea})';
                    else
                        area_aux =[];
                        area_aux(classes) =area_both;
                        valid = valid & (area_aux(spkid)>0)';
                    end
                    spktimes_ = spktimes(valid);
                    spktimes_ = spktimes_-trialtime(itrial,1);

                    spkid_ = spkid(valid);

                    for ineuron=1:length(classes)
                        idx = find(spkid_==classes(ineuron));
                        if isempty(idx)
                            continue
                        end
                        count = histc(spktimes_(idx),timebins);
                        psth(ineuron,itrial,:)=count;
                        psth_smooth(ineuron,itrial,:)=conv(count,rectwin(binwitdh)/binwitdh,'same');
                        if doplot
                            if ~isempty(idx)
                                plot(spktimes_(idx),ineuron,'.k')
                            end
                        end
                    end
                end

                %%
                nlick = nlick_all{ises};
                lick_trace = mean(nlick(:,idx_trace),2);

                right1 = find(trialtype==1 & lick_trace > 0);
                wrong1 = find(trialtype==1 & lick_trace <= 0);
                wrong2 = find(trialtype==2 & lick_trace > 0);
                right2 = find(trialtype==2 & lick_trace <= 0);

                %%

                if computeassemblies==1
                    M = psth_smooth(:,:,periods{isource});
                    M = M(:,:,1:20:end);
                    M = reshape(M,size(M,1),[]);
                    P = assembly_patterns(M);
                else
                    P=patterns_all{iarea,ises,isource};
                end

                psth_proj=zeros(size(P,2),size(psth_smooth,2),length(timebins));
                parfor itrial=1:size(psth_smooth,2)
                    psth_proj(:,itrial,:) = assembly_activity(P,squeeze(psth_smooth(:,itrial,:)));
                end            

                if isempty(P)
                    continue
                end

                % Normalize assembly by PRE activity
                m = squeeze(nanmean(psth_proj(:,:,periods{PRE}),3));
                s = squeeze(nanstd(psth_proj(:,:,periods{PRE}),0,3));
                psth_proj = bsxfun(@rdivide,bsxfun(@minus, psth_proj,m),s);

                mean_assembly_stim_ = squeeze(nanmean(psth_proj(:,:,periods{STIM}),3));
                mean_assembly_trace_ = squeeze(nanmean(psth_proj(:,:,periods{TRACE}),3));
                mean_assembly_rwd_ = squeeze(nanmean(psth_proj(:,:,periods{RWD}),3));

                trialid = [];
                trialid([right1' wrong1'],1)=1;
                trialid([right2' wrong2'],1)=2;
                trialid([right1' right2'],2)=1;
                trialid([wrong1' wrong2'],2)=2;

                mean_assembly_stim_ = mean_assembly_stim_(:,trialorder);
                mean_assembly_trace_ = mean_assembly_trace_(:,trialorder);
                mean_assembly_rwd_ = mean_assembly_rwd_(:,trialorder);
                trialid = trialid(trialorder,:);


                psth_proj = psth_proj(:,trialorder,:);
                nwin=10;
                windows = round(linspace(1, length(trialid),nwin+1));
                mean_assembly_stim1_s =zeros(size(P,2),nwin,length(timebins));
                mean_assembly_stim2_s =zeros(size(P,2),nwin,length(timebins));
                for iwin=1:length(windows)-1
                    mean_assembly_stim1_s(:,iwin,:) = ...
                        squeeze(nanmean(psth_proj(:,intersect(windows(iwin):windows(iwin+1),...
                        find(trialtype(trialorder)==1)),:),2));
                    mean_assembly_stim2_s(:,iwin,:) = ...
                        squeeze(nanmean(psth_proj(:,intersect(windows(iwin):windows(iwin+1),...
                        find(trialtype(trialorder)==2)),:),2));
                end

                %%
                timebins = -1:0.001:1;
                licktime = etime(etype==3);
                if ~issorted(licktime(:,1))
                    licktime = sort(licktime);
                end
                licktime(:,2) = licktime+1;

                %%


                timebins = -1:0.001:1;
                riptime = double(data(ises).RipTS{1})/srate;
                riptime=riptime(:);
                riptime(:,2) = riptime+1;
                rip_prevtrial = zeros(1,size(riptime,1));
                rip_prevtrial_ = zeros(1,size(riptime,1));
                rip_nexttrial = zeros(1,size(riptime,1));
                clear psth psth_smooth
                psth = zeros(length(classes),size(riptime,1),length(timebins));
                psth_smooth = zeros(length(classes),size(riptime,1),length(timebins));

                for irip = 1:size(riptime,1)


                    valid = spktimes>=riptime(irip,1)-1 & spktimes<riptime(irip,2);


                    [v1,idx] = sort(trialtime(:,1));
                    [v2,iidx] = sort(idx);

                    aux = find(riptime(irip,1)<=trialtime(idx(2:end),1)-1 & riptime(irip,1)>trialtime(idx(1:end-1),2)+1);

                    if ~isempty(aux)
                        rip_prevtrial(irip) = find(iidx==aux);

                        rip_prevtrial_(irip) = aux;

                        rip_nexttrial(irip) = find(iidx==aux+1);
                    end


                    if ~isempty(valid_neurons{iarea})
                        area_aux =[];
                        area_aux(classes) =area_;
                        valid = valid & (area_aux(spkid)==valid_neurons{iarea})';
                    else
                        area_aux =[];
                        area_aux(classes) =area_both;
                        valid = valid & (area_aux(spkid)>0)';
                    end
                    spktimes_ = spktimes(valid);
                    spktimes_ = spktimes_-riptime(irip,1);

                    spkid_ = spkid(valid);

                    for ineuron=1:length(classes)
                        idx = find(spkid_==classes(ineuron));
                        if isempty(idx)
                            continue
                        end
                        count = histc(spktimes_(idx),timebins);
                        psth(ineuron,irip,:)=count;
                        psth_smooth(ineuron,irip,:)=conv(count,rectwin(binwitdh),'same');
                        if doplot
    %                         if ~isempty(idx)
    %                             plot(spktimes_(idx),ineuron,'.k')
    %                         end
                        end
                    end
                    %%
                end

                psth_proj_rip=zeros(size(P,2),size(psth_smooth,2),length(timebins));
                parfor irip=1:size(psth_smooth,2)
                    psth_proj_rip(:,irip,:) = assembly_activity(P,squeeze(psth_smooth(:,irip,:)));
                end                        
                %%
                n_rip=zeros(1,length(trialid));
                mean_assembly_rip = zeros(size(P,2),size(trialid,1),length(timebins));
                for irip = unique(rip_prevtrial_)
                    if irip==0
                        continue
                    end
                    idx = find(rip_prevtrial_==irip);
                    mean_assembly_rip(:,irip,:) = nanmean(psth_proj_rip(:,idx,:),2);
                    n_rip(irip) = length(idx);
                end

                % %
                idx = find(timebins>-0.025 & timebins<0.025);
                mean_assembly_rip_ = squeeze(nanmean(mean_assembly_rip(:,:,idx),3));

                nwin=10;
                windows = round(linspace(1, length(trialid),nwin+1));
                mean_assembly_rip_s =zeros(size(P,2),nwin,length(timebins));
                for iwin=1:length(windows)-1
                    mean_assembly_rip_s(:,iwin,:)= nanmean(mean_assembly_rip(:,windows(iwin):windows(iwin+1),:),2);
                end  

                %%
    %             
                %%
                n_stim = {right1 right2 wrong1 wrong2};

                mean_assembly{1,1} = mean_assembly_stim1_s;
                mean_assembly{2,1} = mean_assembly_stim2_s;
    %             mean_assembly{3,1} = mean_assembly_stim;
                mean_assembly{4,1} = mean_assembly_lick_s;
                mean_assembly{5,1} = mean_assembly_rip_s;

                mean_assembly{1,2} = mean_assembly_stim_;
                mean_assembly{2,2} = mean_assembly_trace_;
                mean_assembly{3,2} = mean_assembly_rwd_;
                mean_assembly{4,2} = mean_assembly_lick_;
                mean_assembly{5,2} = mean_assembly_rip_;

                mean_assembly{1,3} = n_stim;
                mean_assembly{2,3} = n_rip;

                mean_assembly_all{iarea,ises,isource} = mean_assembly;

                end

            end
        end
        %%
        save(['Figure5_data2'], 'area_labels',...
            'CA1','PFC','binwitdh','valid_neurons','valid_ses',...
            'mean_assembly_all','seed','-v7.3')

    end

end


%%

timebins = -1:0.001:1;
 
source_label = {'PRE','STIM','TRACE','RWD','ALL'};

CSPLUS=1;
CSMINUS=2;
STIM=3;
LICK=4;
RIP = 5;
R1 = 1;
R2 = 2;
W1 = 3;
W2 = 4;

learn_labels = {'pre-learning' 'post-learning'};

col = [.2 .2 .2;
    .2 .7 .2];

col2 = [.2 .7 .2;
    .6 .9 .6;
    .7 .2 .2;
    .9 .6 .6];

figure(1);clf

for iarea = 1:2
    
    for isource = ALL
    for ilearn = 1:2
        
        %%
        idx = find(learned(:,1)==ilearn-1 & valid_ses(:,iarea));
        ises = idx;

        m_all_rip=[];
        
        for ises_ = ises'

            mean_assembly = mean_assembly_all{iarea,ises_,isource};

            if isempty(mean_assembly)
                continue
            end
            
            m = mean_assembly{RIP,1};
            if size(m,2)==1
                m=m';
            end
            m_all_rip = cat(1,m,m_all_rip);
        end
        
        %%        
        
        norm = m_all_rip(:,:,abs(timebins)>0.05);
        norm =norm(:,:);
        z_m_all_rip = bsxfun(@rdivide, bsxfun(@minus, m_all_rip,nanmean(norm,2)), nanstd(norm,0,2));
        
        idx = find(timebins>-0.025 & timebins<0.025);
        peak = squeeze(nanmean(z_m_all_rip(:,:,idx),3));
        
        sub(4,4,(iarea-1)*2+1,2+(ilearn-1)*2)

        if iarea==1
            edges = -20:0.5:50;
        else 
            edges = -10:0.5:30;
        end
        count = histc(peak(:),edges);
        bar(edges,count,'facecolor','k')
        
        valid = ((nanmean(peak'))>0);
        
        xlabel('Assembly reactivation')
        ylabel('Counts')
        
        title([area_labels{iarea} '-' learn_labels{ilearn} ])
        
        sub(4,4,(iarea-1)*2+(1:2),1+(ilearn-1)*2)
        
        idx = find(timebins>-0.025 & timebins<0.025);
        peak = squeeze(mean(z_m_all_rip(:,:,idx),3));
        peak =peak(find(~isnan(sum(peak,2))),:);
        valid = ((nanmean(peak'))>0);

        aux = nanmean(peak(:,:)')';
        [v, order] = sort(aux(:,end),'descend');
        valid_idx= find(valid);
        
        imagesc(peak((order),:))
        colorbar
        if iarea==1
            caxis([-5 50])
        else
            caxis([-5 20])
        end
        xlabel('Trial block')
        ylabel('Assembly #')
        
        
        sub(4,4,(iarea-1)*2+2,2+(ilearn-1)*2)

        m = nanmean(peak(valid,:));
        s= nanstd(peak(valid,:))./sqrt(sum(~isnan(peak(valid,:))));
        plot(1:10,m,'k','linewidth',2)
        hold on
        plot(1:10,m+s,'k','linewidth',1)
        plot(1:10,m-s,'k','linewidth',1)
        
        m = -nanmean(peak(~valid,:));
        s = nanstd(peak(~valid,:))./sqrt(sum(~isnan(peak(~valid,:))));
        plot(1:10,m,'r','linewidth',2)
        hold on
        plot(1:10,m+s,'r','linewidth',1)
        plot(1:10,m-s,'r','linewidth',1)
        
        if iarea==1
            ylim([-.5 20])
        else
            ylim([-.5 10])
        end
        xlim([.5 11.5])
        
        ylabel('Mean assembly reactivation (+SEM)')
        xlabel('Trial block')
        
    end
    end
end


function avg_deconPlot_cc(nSIDs,mode)

% average_deconfiles;


% load(fullfile(datafolder,'avg_decon.mat'));

% OR

load(fullfile(datafolder,'avg_decon_deconv.mat'));

rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
%% Plot CC
for ri = 1%:length(rois)
    lcc = squeeze(avg_decon.cc(ri,:,:,:));
    lcc_s = squeeze(avg_decon.ccs(ri,:,:,:));
    conidx = avg_decon.conidx;
    cohidx = avg_decon.cohidx;
    decon.ehdr = squeeze(lcc(1,:,:));
    decon.ehdrste = squeeze(lcc_s(1,:,:));
    model.ehdr = squeeze(lcc(2,:,:));
    decon.time=0.25:0.5:40.5;
    model.time=decon.time;
    clist = brewermap(10,'PuOr');
    f= figure; hold on
    % plot #1: contrast against 0% coherence (no change)
    lconidx = find(cohidx==0);
    convalues = [0.25 0.5 0.75 1];
    colmap = [4 3 2 1];
%     subplot(2,1,1), hold on
    if isempty(strfind(mode,'nomodel'))
        for i = 1:length(convalues)
            ci = colmap(find(convalues==conidx(lconidx(i)),1));
            plot(model.time,model.ehdr(lconidx(i),:),'Color',clist(ci,:));
        end
    end
    if strfind(mode,'nozero')
        for i = 2:length(convalues)
            ci = colmap(find(convalues==conidx(lconidx(i)),1));
            errbar(decon.time,decon.ehdr(lconidx(i),:),decon.ehdrste(lconidx(i),:),'Color',clist(ci,:));
            p(i-1) = plot(decon.time,decon.ehdr(lconidx(i),:),'o','MarkerSize',3,'MarkerFaceColor',clist(ci,:),'MarkerEdgeColor',[1 1 1]);
        end
        legs = {'\Delta Con +25%','\Delta Con +50%','\Delta Con +75%'};
    else
        for i = 1:length(convalues)
            ci = colmap(find(convalues==conidx(lconidx(i)),1));
            plot(decon.time,decon.ehdr(lconidx(i),:),'o','MarkerSize',3,'MarkerFaceColor',clist(ci,:),'MarkerEdgeColor',[1 1 1]);
            errbar(decon.time,decon.ehdr(lconidx(i),:),decon.ehdrste(lconidx(i),:),'Color',clist(ci,:));
        end
        legs = {'Contrast=+0%','Contrast=+25%','Contrast=+50%','Contrast=+75%'};
    end
    if ri==8
        l = legend(p,legs,'FontSize',7,'FontName','Helvetica');
        set(l,'box','off');
%         v = get(l,'title');
%         set(v,'string','\Delta Con');
    end
    % plot #2: coherence against 25% contrast (no change)
    lcohidx = find(conidx==0.25);
    a = axis;
    if ri==1
        axis([0 15 -0.75 3.5]);
    else
        axis([0 15 -0.5 1.5]);
    end
    set(gca,'XTick',[0 5 10 15]);
    set(gca,'YTick',[-0.5 0 1]);
    xlabel('Time (s)');
    ylabel('\Delta signal (%)');
    drawPublishAxis('figSize=[4.45,3]');
    
%     if any(ri==[1 8])
%     savepdf(f,sprintf('~/proj/att_awe/talks/data_figures/avg_cc_con_%s.pdf',rois{ri}));
%     end
    savepdf(f,fullfile(datafolder,'avg_decon',sprintf('avg_cc_con_%s',rois{ri})));
    f = figure; hold on
%     subplot(2,1,2), hold on
    if isempty(strfind(mode,'nomodel'))
        if ~isempty(strfind(mode,'nozero'))
            for i = 2:5
                plot(model.time,model.ehdr(lcohidx(i),:),'Color',clist(i+5,:));
            end
        else
            for i = 1:5
                plot(model.time,model.ehdr(lcohidx(i),:),'Color',clist(i+5,:));
            end
        end
    end
    
    if ~isempty(strfind(mode,'nozero'))
        for i=2:5            
            errbar(decon.time,decon.ehdr(lcohidx(i),:),decon.ehdrste(lcohidx(i),:),'Color',clist(i+5,:));

            p(i-1) = plot(decon.time,decon.ehdr(lcohidx(i),:),'o','MarkerSize',3,'MarkerFaceColor',clist(i+5,:),'MarkerEdgeColor',[1 1 1]);
        end
        legs = {'\Delta Coh +25%','\Delta Coh +50%','\Delta Coh +75%','\Delta Coh +100%'};
    else
        for i=1:5
            plot(decon.time,decon.ehdr(lcohidx(i),:),'o','MarkerSize',3,'MarkerFaceColor',clist(i+5,:),'MarkerEdgeColor',[1 1 1]);
            errbar(decon.time,decon.ehdr(lcohidx(i),:),decon.ehdrste(lcohidx(i),:),'Color',clist(i+5,:));
        end
        legs = {'Coherence+0%','Coherence+25%','Coherence+50%','Coherence+75%','Coherence+100%'};
    end
    if ri==8
        l = legend(p,legs,'FontSize',7,'FontName','Helvetica');
        set(l,'box','off');
    end
    
    a = axis;
    axis([0 15 -0.5 1.5]);
    set(gca,'XTick',[0 5 10 15]);
    set(gca,'YTick',[-0.5 0 1]);
    
    xlabel('Time (s)');
    ylabel('\Delta signal (%)');

    drawPublishAxis('figSize=[4.45,3]');
    
%     if any(ri==[1 8])
%     savepdf(f,sprintf('~/proj/att_awe/talks/data_figures/avg_cc_coh_%s.pdf',rois{ri}));
%     end
    savepdf(f,fullfile(datafolder,'avg_decon_allf',sprintf('avg_cc_coh_%s',rois{ri})));
end

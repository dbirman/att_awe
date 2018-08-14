function avg_deconPlot_time(nSIDs,mode)

% average_deconfiles;


% load(fullfile(datafolder,'avg_decon.mat'));

% OR

load(fullfile(datafolder,'avg_decon_deconv.mat'));

rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};


%% plot time (only 100% 100%)

% max = [3 1.5];
% min = [-1 -0.5];
for ri = 1%:length(rois)
    % setup
    ltime = squeeze(avg_decon.time(ri,:,:,:));
    ltime_s = squeeze(avg_decon.times(ri,:,:,:));
%     cc_s = squeeze(cc_s(ri,:,:,:));
    conidxs = avg_decon.timconidx;
    cohidxs = avg_decon.timcohidx;
    timidxs = avg_decon.timidxs;
    decon.ehdr = squeeze(ltime(1,:,:));
    decon.ehdrste = squeeze(ltime_s(1,:,:));
    model.ehdr = squeeze(ltime(2,:,:));
    decon.time=0.25:0.5:40.5;
    model.time=decon.time;
    % plotting
    clist = brewermap(6,'Greys');
    clist(1:5,:) = clist(2:6,:);
    f = figure; hold on
    % first plot will be contrast timing, when contrast goes up to 50% or 100%,
    % we'll draw each of these with increasing contrast colors
    conopts = [0.50 0.50 1 1];
    cohopts = [0.25 1 0.25 1];
    flip = [0.5 1 2 4 8];
    colpos = [1 2 3 4 5];
    for sub = 4
%         subplot(2,2,sub); hold on
%         title(sprintf('Con: %i%% Coh: %i%%',conopts(sub)*100,cohopts(sub)*100));
        idxs = find(logical(conidxs==conopts(sub)).*logical(cohidxs==cohopts(sub)));
        for ii = 1:length(idxs)
            i = idxs(ii);
            ci = colpos(find(flip==timidxs(i),1));
            errbar(decon.time,decon.ehdr(i,:),decon.ehdrste(i,:),'Color',clist(ci,:));
            p(ii) = plot(decon.time,decon.ehdr(i,:),'o','MarkerSize',3,'MarkerFaceColor',clist(ci,:),'MarkerEdgeColor',[1 1 1]);
            if isempty(strfind(mode,'nomodel'))
                plot(model.time,model.ehdr(i,:),'Color',clist(ci,:));
            end
        end
        if ri==8
%             legend(p,{'250 ms','500 ms','1000 ms','2000 ms','4000 ms'},'FontSize',7,'FontName','Helvetica');
        end
%         xlabel('Time (s)');
%         ylabel('\Delta signal (%)');
        a = axis;
        if ri==1
            axis([0 15 -0.75 3.5]);
        else
            axis([0 15 -0.5 1.5]);
        end
        set(gca,'XTick',[0 5 10 15],'YTick',[-0.5 0 1]);
    end
    drawPublishAxis('figSize=[3.5,3.5]');
%     if any(ri==[1 8])
%     savepdf(f,sprintf('~/proj/att_awe/talks/data_figures/avg_time_100100_%s.pdf',rois{ri}));
%     end
    savepdf(f,fullfile(datafolder,'avg_decon',sprintf('avg_time_100100_%s',rois{ri})));
end

%% plot time

% These plots are great for visualizing, but not great for a paper. Just
% plot one example (100% 100% in a larger plot instead)
% max = [3 1.5];
% min = [-1 -0.5];
% for ri = 1:length(rois)
%     % setup
%     ltime = squeeze(avg_decon.time(ri,:,:,:));
%     ltime_s = squeeze(avg_decon.times(ri,:,:,:));
% %     cc_s = squeeze(cc_s(ri,:,:,:));
%     conidxs = avg_decon.timconidx;
%     cohidxs = avg_decon.timcohidx;
%     timidxs = avg_decon.timidxs;
%     decon.ehdr = squeeze(ltime(1,:,:));
%     decon.ehdrste = squeeze(ltime_s(1,:,:));
%     model.ehdr = squeeze(ltime(2,:,:));
%     decon.time=0.25:0.5:40.5;
%     model.time=decon.time;
%     % plotting
%     clist = brewermap(6,'Greys');
%     clist(1:5,:) = clist(2:6,:);
%     f = figure; hold on
%     % first plot will be contrast timing, when contrast goes up to 50% or 100%,
%     % we'll draw each of these with increasing contrast colors
%     conopts = [0.50 0.50 1 1];
%     cohopts = [0.25 1 0.25 1];
%     flip = [0.5 1 2 4 8];
%     colpos = [1 2 3 4 5];
%     for sub = 1:4
%         subplot(2,2,sub); hold on
%         title(sprintf('Con: %i%% Coh: %i%%',conopts(sub)*100,cohopts(sub)*100));
%         idxs = find(logical(conidxs==conopts(sub)).*logical(cohidxs==cohopts(sub)));
%         for i = idxs
%             ci = colpos(find(flip==timidxs(i),1));
%             plot(decon.time,decon.ehdr(i,:),'o','MarkerSize',7,'MarkerFaceColor',clist(ci,:),'MarkerEdgeColor',[1 1 1]);
%             errbar(decon.time,decon.ehdr(i,:),decon.ehdrste(i,:),'Color',clist(ci,:));
%             if isempty(strfind(mode,'nomodel'))
%                 plot(model.time,model.ehdr(i,:),'Color',clist(ci,:));
%             end
%         end
%         xlabel('Time (s)');
%         ylabel('Response (%signal/s)');
%         a = axis;
%         axis([0 15 a(3) a(4)]);
%     end
%     savepdf(f,fullfile(datafolder,'avg_decon',sprintf('avg_time_%s',rois{ri})),'subplots=[2,2]','figsize=[8,8]');
% end
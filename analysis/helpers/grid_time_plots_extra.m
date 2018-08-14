%% Grid plot
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
x = 0.25:0.5:40.5;
riopts = [2:7];
for ri = 1:length(riopts)
    ropt = riopts(ri);
    % roi | model | index | value
    load(fullfile(datafolder,'avg_decon.mat'));

    % get V1, real data, and squeeze
    cc = squeeze(avg_decon.cc(ropt,1,:,:));
    ccs = squeeze(avg_decon.ccs(ropt,1,:,:));

    conidx = avg_decon.conidx;
    cohidx = avg_decon.cohidx;
    ucon = unique(conidx);
    ucoh = unique(cohidx);
    h = figure;
    cmap = brewermap(11,'PuOr');
    colors = repmat([0,0,0],20,1);
    colors(17,:) = cmap(7,:);
    colors(18,:) = cmap(8,:);
    colors(19,:) = cmap(9,:);
    colors(20,:) = cmap(10,:);
    colors(11,:) = cmap(4,:);
    colors(6,:) = cmap(3,:);
    colors(1,:) = cmap(2,:);
    cons = {'+0%','+25%','+50%','+75%'};
    cohs = cons; cohs{end+1} = '+100%';
    for conii = 1:length(ucon)
        coni = 5-conii;
        for cohi = 1:length(ucoh)
            idx = logical((conidx==ucon(coni)).*(cohidx==ucoh(cohi)));
            i = (conii-1)*length(ucoh)+cohi;
            subplot(4,5,i)
    %         errbar(x,cc(idx,:),ccs(idx,:));
%             plot(x,cc(idx,:),'-','Color',colors(i,:),'LineWidth',1);
            plot(x,cc(idx,:),'-k','LineWidth',1);
            if ri==1
                axis([0 20 -1 3]);
            else
                axis([0 20 -1 2]);
            end
            if ropt==2
                set(gca,'XTick',[0 10]);
                set(gca,'YTick',[0 1]);
            else
                set(gca,'XTick',[0]);
                set(gca,'YTick',[0]);
            end
%             if coni==2 && cohi==2
%                 xlabel('Time (s)');
%                 ylabel('\Delta signal (%)');
%             end
%             if cohi==1
%                 ylabel(cons{coni});
%             end
%             if coni==1
%                 xlabel(cohs{cohi});
%             end
            set(gca,'Visible','off');
            drawPublishAxis('figSize=[4.5,4.5]');
        end
    end
%     savepdf(h,fullfile('~/proj/att_awe/talks/data_figures/',sprintf('grid_plot_%s.pdf',rois{ropt})));
    savepdf(h,fullfile(datafolder,'grid_extra',sprintf('grid_plot_%s.pdf',rois{ropt})));
    %%%%%%%%%%%%%%% OLD CODE %%%%%%%%%%%%%%%%%

%     %% Grid time plot
%     x = 0.25:0.5:40.5;
% 
%     h = figure;
%     time = squeeze(avg_decon.time(ropt,1,:,:));
%     gdeltas = {'+50/+25','+50/+100','+100/+25','+100/+100'};
%     lengths = {'250 ms','500 ms','1000 ms','2000 ms','4000 ms'};
%     for group = 1:4
%         for ti = 1:5
%             idx = (group-1)*5+ti;
%             subplot(4,5,idx);
%             plot(x,time(idx,:),'-k','LineWidth',1);
%             if ri==1
%                 axis([0 20 -1 3]);
% %                 set(gca,'YTick',[0 2 4],'YTickLabel',[0 2 4]);
%             else
%                 axis([0 20 -1 2]);
%             end
%             if ti==1
%                 ylabel(gdeltas{group});
%             end
%             if group==4
%                 xlabel(lengths{ti});
%             end
%             set(gca,'XTick',[0 10],'YTick',[0 1]);
%             drawPublishAxis('figSize=[4.5,4.5]');
%         end
%     end
% %     savepdf(h,fullfile('~/proj/att_awe/talks/data_figures/',sprintf('grid_timeplot_%s.pdf',rois{ropt})));
%     savepdf(h,fullfile(datafolder,'grid_extra',sprintf('grid_timeplot_%s.pdf',rois{ropt})));
end

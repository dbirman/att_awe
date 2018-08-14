%% 3-D grid plot (For Figure 1)
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
load(fullfile(datafolder,'avg_decon.mat'));

cmap = brewermap(11,'PuOr');
riopts = 1:8;

for ri = 1:length(riopts)
    ropt = riopts(ri);
    h = figure; hold on
    % pull just V1
    cc = squeeze(avg_decon.cc(ropt,1,:,:));
    time = squeeze(avg_decon.time(ropt,1,:,:));

    x = 0.25:0.5:15.5;
    ucon = [0.25 0.5 0.75 1];
    ucoh = [0 0.25 0.5 0.75 1];
    utime = [0.5 1 2 4 8];
    xtoff = [-0.8 -0.6 -0.4 -0.2 0.2]*max(x)/2;
    ytoff = [-0.8 -0.6 -0.4 -0.2 0.2]/5;% + [-0.05 -0.05 -0.05 -0.05 .05];

    % now plot the time data
    clist = brewermap(6,'Greys');
    clist(1:5,:) = clist(2:6,:);
    for cond = 1:20
        ccon = avg_decon.timconidx(cond);
        ccoh = avg_decon.timcohidx(cond);
        ctim = avg_decon.timidxs(cond);
        yidx = ((find(ucon==ccon,1)-1)*3);
        xidx = (find(ucoh==ccoh,1)-1)*(max(x)+5);
        p = plot(xtoff(find(utime==ctim,1))+xidx+x,ytoff(find(utime==ctim,1))+(yidx+time(cond,1:length(x)))/10,'-','Color',clist(find(utime==ctim,1),:));

    end
    
    % first plot the cohcon data at time = 2.5
    for cond = 2:20
        ccon = avg_decon.conidx(cond);
        ccoh = avg_decon.cohidx(cond);
        yidx = ((find(ucon==ccon,1)-1)*3);
        xidx = (find(ucoh==ccoh,1)-1)*(max(x)+5);
        p = plot(xidx+x,(yidx+cc(cond,1:length(x)))/10,'-k');
    end
    
    plot([0 0],[0 0.1]-.25,'-k');
    plot([0 10],[-0.25 -0.25],'-k');
    
    set(gca,'XColor',[1 1 1]);
    set(gca,'YColor',[1 1 1]);
    drawPublishAxis('figSize=[5.25,5]');
    savepdf(h,fullfile(datafolder,'grid',sprintf('grid_plot_%s.pdf',rois{ropt})));
end

% set(gca,'XTick',
% % add helper lines
% 
% % put a line at 2500 ms
% % plot3(1:4
% 
% set(gca,'XTick',[0.5 1.5 2.5 3.5]*max(x)+[0 5 10 15],'XTickLabel',ucon*100);
% xlabel('Contrast change (%)');
% 
% set(gca,'YTick',[1 2 3 4 5],'YTickLabel',utime/2*1000);
% ylabel('Time (ms)');
% 
% set(gca,'ZTick',[1 2 3 4 5]*4/10,'ZTickLabel',ucoh*100);
% zlabel('Coherence change (%)');

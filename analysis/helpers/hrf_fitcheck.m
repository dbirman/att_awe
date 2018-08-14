%% HRF Fit check

% for each subject, try fitting an HRF individually to every ROI

% s0300
subj = 's0300';
folders = {'s030020160404','s030020160411','s030020160509','s030020160528','s030020160529'};
datas = {};
for fi = 1:length(folders)
    load(fullfile(sprintf('C:/Users/dan/proj/COHCON_DATA/%s_data.mat',folders{fi})));
    datas{end+1}=data;
end
load(fullfile(sprintf('C:/Users/Dan/proj/COHCON_DATA/%s_fithrf.mat',subj)));
fits = cell(size(fithrf.ROIs));
for ri = 1:length(fithrf.ROIs)
    disp(sprintf('Starting %s',fithrf.ROIs{ri}));
    tdatas = datas;
    for di = 1:length(datas)
        tdatas{di}.rtSeries = tdatas{di}.rtSeries(ri);
        tdatas{di}.tSeries = tdatas{di}.tSeries(ri);
        tdatas{di}.ROIs = tdatas{di}.ROIs(ri);
    end
    fits{ri} = fitCCTimecourseROIModel(tdatas,'fithrf');
    clear tdatas
    disp(sprintf('Finished %s',fithrf.ROIs{ri}));
end

%% plot against the canonical HRF fit for all ROIs
x = 0.25:.5:49.75;
canon = cc_doublegamma(x,fithrf.params);
hrfs = zeros(length(fits),length(x));
for hi = 1:size(hrfs,1)
    cur = cc_doublegamma(x,fits{hi}.params);
    cur = cur ./ max(cur);
    hrfs(hi,:) = cur;
end

canon = canon ./ max(canon);
%% plot
figure, hold on
plot(x,hrfs);
plot(x,canon,'-k','LineWidth',4);
title('Comparison of Normalized Individual vs. Across ROI HRF Fits');
xlabel('Time');
ylabel('Amplitude (a.u.)');
drawPublishAxis


fname = fullfile('C:/Users/Dan/proj/COHCON_DATA','avg_hrf.pdf');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'InvertHardCopy','off');
set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fname,'-dpdf');
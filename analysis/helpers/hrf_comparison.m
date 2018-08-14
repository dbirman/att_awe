%% Load (in hrfs)

load(fullfile(datafolder,'avg_hrf.mat'));
    
%% plot hrfs for comparison
figure;
for ri = 1:4
    subplot(4,1,ri); hold on
    plot(squeeze(hrfs(:,ri,:))');
    title(strrep(dataopts{ri},'_','-'));
end

%% Plot only HRFs from the top-25 condition
h = figure; hold on

t = 0.25:.5:40.5;
for si = 1:size(hrfs,1)
    cdat = squeeze(hrfs(si,4,:));
    cdat = cdat./max(cdat(:));
    plot(t,cdat,'-k','LineWidth',1);
end
avg = squeeze(mean(hrfs(:,4,:)));
avg = avg ./ max(avg);
plot(t,avg,'--r','LineWidth',1);
axis([0 40 -0.5 1.25]);
xlabel('Time (s)');
ylabel('Signal change (%)');

drawPublishAxis('figSize=[8.9,5]');

savepdf(h,fullfile(datafolder,'avg_fmri','avg_hrfs.pdf'));
%% compare top-2 and top-25 HRFs
figure
plot((squeeze(hrfs(:,3,:))-squeeze(hrfs(:,4,:)))');

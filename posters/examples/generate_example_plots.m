%% Plot the neural response example

f = figure;
hold on

x = 0:.01:1;
p.conmodel = 3;
p.conalpha = 10;
p.conkappa = 2;
y = conModel(x,p);

cmap = brewermap(7,'PuOr');

plot(x*100,y,'Color',cmap(2,:));
plot(x*100,y*0.5,'Color',cmap(6,:));

legend('Contrast','Coherence','Location','northwest');

xlabel('Stimulus Strength (%)');
ylabel('Neural Firing Rate (a.u.)');

set(gca,'XTick',[0 100]);
set(gca,'YTick',[0 9],'YTickLabel',{'',''});

drawPublishAxis
savepdf(f,fullfile('~/proj/att_awe/posters/examples/neural_response.pdf'));

%% Plot HRF response example

f = figure;

y = spm_hrf(0.5);
y = y/max(y);

plot(y,'-k');

xlabel('Time (s)');
ylabel('BOLD Signal (%)');

set(gca,'XTick',[0 65],'XTickLabel',{'',''});
set(gca,'YTick',[0 1],'YTickLabel',{'',''});
drawPublishAxis
savepdf(f,fullfile('~/proj/att_awe/posters/examples/hrf_response.pdf'));

%% Timeseries example

events = [zeros(1,4) 1 zeros(1,5) 0.75 zeros(1,5) 0.5 zeros(1,10) 1.5 zeros(1,5)];
hrf = spm_hrf(0.5);
y = conv(events,hrf)';
y = y + randn(1,length(y))/100;
y = y/max(y);

f = figure;



plot(y,'-k');
xlabel('Time (s)');
ylabel('BOLD Signal (%)');
set(gca,'XTick',[0 length(y)],'XTickLabel',{'',''});
set(gca,'YTick',[0 1],'YTickLabel',{'',''});
drawPublishAxis

savepdf(f,fullfile('~/proj/att_awe/posters/examples/timeseries.pdf'));

%%
h = figure;
x = -3:.01:3;
y = normpdf(x);
plot(x,y);
drawPublishAxis
savepdf(h,fullfile('~/proj/att_awe/posters/examples/normal.pdf'));
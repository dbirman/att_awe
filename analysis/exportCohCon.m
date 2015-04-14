%% To start off, let's get files

mglSetSID('s300')

global analysis

[analysis.datFolder, analysis.anFolder] = getSubjDataFolder('cohcon',false);

year = date;
year = year(end-1:end);

files = dir(sprintf('%s/%s*mat',analysis.datFolder,year));

%% Loop over files, fixing them for analysis

fixFiles(files);

%% Loop over files and load them into expHolder

expHolder = loadExp(files);

%% Loop over exp files and export them to csv

for ei = 1:length(expHolder)
%     try
        cohCon2csv(expHolder{ei},true);
%     catch
%         disp(sprintf('Experiment file %i not generated...',ei));
%     end
end

%% The following functions use the most recent stimulus file

load(fullfile(analysis.datFolder,files(end).name));

%% Generate staircase graphs

% staircaseplots(stimulus);

%% Send main to CSV file
plotting = cohCon_discFuncs(stimulus,1);
cohCon_plo2csv(plotting);

%% Send catch to CSV file
cat = cohCon_catPerf(stimulus,0);
per2csv(cat);

%%
high = out.fit.x(find(out.fit.y>.99,1,'first'));
if isempty(high)
    x = .001:.001:.999;
    y = weibull(x,out.fit.fitparams);
    high = x(find(y>.99,1,'first'));
end
if isempty(high)
    high = x(find(y==y(end),1,'first'));
end
high
low = out.fit.x(find(out.fit.y<.51,1,'last'));
if isempty(low)
    x = .001:.001:.999;
    y = weibull(x,out.fit.fitparams);
    low = x(find(y<.51,1,'last'));
end
if isempty(low)
    low = x(find(y==y(1),1,'last'));
end
low

%% Plotting a figure with weibull fits across a range of values

f = figure
hold on
colors = {'-r','-g','-b','-c'}; vals = {};
for i = 1:4
    figure(f)
    out = doStaircase('threshold',stimulus.staircase{2,i},'type','weibull');
    plot(out.fit.x,out.fit.y,colors{i},'LineWidth',3);
    vals{i} = num2str(stimulus.pedestals.contrast(i));
end
legend(vals);
axis([0 .15 .475 1.025]);
set(gca,'fontsize',24);
set(gca,'LineWidth',3);

%% Plotting the %correct for the weibull fit

figure
set(gcf,'renderer','painters')
hold on
plot(out.fit.signal,out.fit.pcorrect,'*b');
plot(out.fit.x,out.fit.y)
a = axis();
axis([a(1) a(2) .48 1.02]);
set(gca,'fontsize',24);
set(gca,'LineWidth',3);
set(gca,'XTick',[0 .1 .2]);
set(gca,'YTick',0:.1:1);

%% Plotting both attention condition's weibull fits

for i = 1:1
    figure
    hold on    
    out = doStaircase('threshold',stimulus.staircase{2,i},'type','weibull');
    outc = doStaircase('threshold',stimulus.stairCatch{2,i},'type','weibull');
    
    plot(out.fit.signal,out.fit.pcorrect,'ob','LineWidth',3,'MarkerFaceColor','b');
    plot(out.fit.x,out.fit.y,'-b');
    
    plot(outc.fit.signal,outc.fit.pcorrect,'or','LineWidth',3,'MarkerFaceColor','r');
    plot(outc.fit.x,outc.fit.y,'-r');
    hold off
    axis([0 .2 .48 1.02]);
set(gca,'fontsize',24);
set(gca,'LineWidth',3);
end

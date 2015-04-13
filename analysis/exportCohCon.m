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

%%

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
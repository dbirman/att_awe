%% To start off, let's get files

mglSetSID('s300')

global analysis

[analysis.datFolder, analysis.anFolder] = getSubjDataFolder;

year = date;
year = year(end-1:end);

files = dir(sprintf('%s/%s*mat',analysis.datFolder,year));

%% Loop over files, fixing them for analysis

fixFiles(files);

%% Loop over files and load them into expHolder

expHolder = loadExp(files);

%% Loop over exp files and export them to csv

for ei = 1:length(expHolder)
    try
        exp2csv(expHolder{ei});
    catch
        disp(sprintf('Experiment file %i not generated...',ei));
    end
end

%% The following functions use the most recent stimulus file

load(fullfile(analysis.datFolder,files(end).name));

%% Generate staircase graphs

staircaseplots(stimulus);

%% Generate discrimination functions

% The idea here is just to compare the single and dual task performance on
% the same graphs.
plotting = gen_discFuncs(stimulus,0);

%% Send plotting to CSV file

plo2csv(plotting);

% % % % % % % %% normalize plotting
% % % % % % % 
% % % % % % % [regplotting, rplottingsd] = normPlotting(plotting,0);
% % % % % % % [normplotting, nplottingsd] = normPlotting(plotting,1);
% % % % % % % 
%% plot disc
% gen_discPlots(plotting,stimulus,0);
% % % % % % % gen_discPlots(regplotting,rplottingsd,stimulus);
% % % % % % % gen_discPlots(normplotting,nplottingsd,stimulus);

% % % % % % % % %% Generate peripheral performance
% % % % % % % % 
% % % % % % % % peripheral = gen_perPerf(stimulus);

%% Send peripheral to CSV file

per2csv(peripheral);

% % % % % % % % %% Plot peripheral
% % % % % % % % 
% % % % % % % % gen_perPlots(peripheral,0);
% % % % % % % % gen_perPlots(peripheral,1);
% % % % % % % % 
% % % % % % % % %% Generate Performance Plots
% % % % % % % % 
% % % % % % % % % The idea here is to have a plot that shows the dual task performance in
% % % % % % % % % comparison with the single task performance for both gender and
% % % % % % % % % contrast/noise at the same time.
% % % % % % % % 
% % % % % % % % gen_perf(stimulus,normplotting,peripheral);
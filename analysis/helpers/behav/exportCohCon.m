%% To start off, let's get files
function allData = exportCohCon(allData)

doEye = false;
scan = false;

global analysis

[analysis.datFolder, analysis.anFolder] = getSubjDataFolder('cohcon',scan);

year = date;
year = year(end-1:end);

files = dir(sprintf('%s/%s*mat',analysis.datFolder,year));

%% Loop over files, fixing them for analysis

fixFiles(files);

%% Loop over files and load them into expHolder

expHolder = loadExp(files);

%% Loop over exp files and export them to csv

for ei = 1:length(expHolder)
    %     try                        %EYE SKIP
    cohCon2csv(expHolder{ei},~doEye);
    %     catch
    %         disp(sprintf('Experiment file %i not generated...',ei));
    %     end
end
if doEye
    return
end

%% The following functions use the most recent stimulus file

load(fullfile(analysis.datFolder,files(end).name));

%% If Scan: Concatenate across stimulus runs
if ~isempty(strfind(analysis.anFolder,'scan'))
    taskCued = {};
    taskMiscued = {};
    control = {};
    for i = 1:length(expHolder)
        expi = expHolder{i}{1};
        for j = 1:2
            taskCued{i,j} = expi.stimulus.staircase{j};
            taskMiscued{i,j} = expi.stimulus.stairCatch{j};
            for k = 1:size(expi.stimulus.nocatchs.staircase,2)
                
%                 try
                control{i,j,k} = expi.stimulus.nocatchs.staircase{j,k};
%                 catch
%                 end
            end
        end
    end
end

%% Build up dataset for Performance vs. Stimulus Intensity

I = struct;

I.catch.con = [];
I.nocatch.con = [];
I.catch.coh = [];
I.nocatch.coh = [];
I.catch.side = [];
I.nocatch.side = [];
I.catch.task = [];
I.nocatch.task = [];
I.catch.catch = [];

% What we want to do is get every single trial's left-right difference for
% contrast and coherence, and which side was chosen. Then we can bin and
% figure out
for i = 1:length(expHolder)
    expi = expHolder{i}{1};
    if length(unique(expi.randVars.task))==2
        % catch
        I.catch.coh = [I.catch.coh expi.randVars.rCoh-expi.randVars.lCoh];
        I.catch.con = [I.catch.con expi.randVars.rCon-expi.randVars.lCon];
        I.catch.side = [I.catch.side expi.response];
        I.catch.task = [I.catch.task expi.randVars.task];
        I.catch.catch = [I.catch.catch expi.parameter.catch];
    else
        % nocatch
        I.nocatch.coh = [I.nocatch.coh expi.randVars.rCoh-expi.randVars.lCoh];
        I.nocatch.con = [I.nocatch.con expi.randVars.rCon-expi.randVars.lCon];
        I.nocatch.side = [I.nocatch.side expi.response];
        I.nocatch.task = [I.nocatch.task expi.randVars.task];
    end
end


%% Send main to CSV file
[plotting, mfits, maf] = cohCon_discFuncs(stimulus,0);
cohCon_plo2csv(plotting);

%% Send catch to CSV file
[cat, cfits, caf] = cohCon_catPerf(stimulus,0);
per2csv(cat);

%% Send nocatch to CSV file
[nocat, nfits, nocaf] = cohCon_nocatPerf(stimulus,0);
nocat2csv(nocat);

%% Send fits to CSV file
cohCon_fits2csv(nfits,mfits,cfits);

%% Eye data

eye = {};
for ei = 1:length(expHolder)
    eye{ei} = expHolder{ei}{2};
end

%% Save into allData

allData.behav.I = I;
allData.behav.maf = maf;
allData.behav.caf = caf;
allData.behav.nocaf = nocaf;
allData.behav.mfits = mfits;
allData.behav.nfits = nfits;
allData.behav.cfits = cfits;
allData.behav.eye = eye;
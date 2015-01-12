%% To start off, let's get files

datFolder = getSubjDataFolder;

year = date;
year = year(end-1:end);

files = dir(sprintf('%s/%s*mat',datFolder,year));

%% Loop over files, fixing them for analysis

for fi = 1:length(files)
    curFile = sprintf('%s/%s',datFolder,files(fi).name);
    fixFile = sprintf('%s/%s',datFolder,['f' files(fi).name]);
    if isfile(fixFile)
        disp(sprintf('(noisecon_analysis) Found fixed file for %s, ignoring...',curFile));
    else
        success = genExp(curFile,fixFile);
        if success
            disp(sprintf('(noisecon_analysis) Output in: %s',fixFile));
        else
            warning(sprintf('(noisecon_analysis) File %s was not fixed!!!!',curFile));
        end
    end
end

%% Loop over files and load them into expHolder

expHolder = {};
for fi = 1:length(files)
    fixFile = sprintf('%s/%s',datFolder,['f' files(fi).name]);
    load(fixFile);
    expHolder{end+1} = exp;
end

%% ANALYSIS TIME

% LIST OF ANALYSES THAT I EXPECT TO BE IMPORTANT:
%   -   Do errors during a dual trial cause the trial to fail
%       a) Are errors during the first half vs second half more important?
%   -   Is reaction time different?
%       a) Male vs. Female vs. None peripheral

%% Check for reaction time differences

%               type = M1 F2 O3
% using gender(single/dual,trial,type,RT)

for i = 1:length(expHolder)
    cExp = expHolder{i};
    main = exp{1};
    per = exp{2};
    rVars = exp{3}.runVars;
    
    % Peripheral Task Differences
    if rVars.dual(i) == 1
        % This is a dual run that is loaded
        for t = 1:length(per.trials)
            if per.randVars.respond==0
                type = 3;
            else
                type = per.
            gender(2,t,
        end
    else
        % This is a non-dual run
    end
end

%% Check for errors during dual trials

for i = 1:length(expHolder)
    cExp = expHolder{i};
    main = exp{1};
    per = exp{2};
    rVars = exp{3}.runVars;
    if rVars.dual(i) == 1
        % This is a dual run that is loaded
        
    end
end

%% Generate staircase graphs

% Depending on the type of the loaded run, generate the staircase figures
stairtype = 'dualstaircase';

dispTexts = {'Noise','Contrast'};
for num = 1:2
    plotting = zeros(2,3);
    
    dispText = dispTexts{num};
    
    try
        figure % this is the 'staircase' figure
        title(sprintf('%s, DUAL Staircase plot (R->G->B high)',dispText));
        hold on
        drawing = {'-r' '-g' '-b'
            '--r' '--g' '--b'};
        for cues = 1:2
            for ped = 1:3
                testV = [];
                % for each cue/ped combo, this pulls out and adds to the test
                % value list. This will be used to display the 'staircases'.
                for i = 1:length(stimulus.(stairtype){num,cues,ped})
                    testV = [testV stimulus.(stairtype){num,cues,ped}(i).testValues];
                end
                plot(testV,drawing{cues,ped});
            end
        end
    catch
        warning('Dual task staircases failed');
    end
    stairtype = 'staircase';
    try
        figure % this is the 'staircase' figure
        title(sprintf('%s, Staircase plot (R->G->B high)',dispText));
        hold on
        drawing = {'-r' '-g' '-b'
            '--r' '--g' '--b'};
        for cues = 1:2
            for ped = 1:3
                testV = [];
                % for each cue/ped combo, this pulls out and adds to the test
                % value list. This will be used to display the 'staircases'.
                for i = 1:length(stimulus.(stairtype){num,cues,ped})
                    testV = [testV stimulus.(stairtype){num,cues,ped}(i).testValues];
                end
                plot(testV,drawing{cues,ped});
            end
        end
        if stimulus.p.staircase(end).trialNum == 0
            % Ignore the last staircase, it was just reset and has no trials
            % (and causes an error)
            usestair = stimulus.p.staircase(1:end-1);
        else
            usestair = stimulus.p.staircase;
        end
    catch
        warning('Single task staircases failed');
    end
end
% Peripheral task
try
    disp('Gender performance during NOISE task');
    doStaircase('threshold',stimulus.p.dualstaircase{1},'dispFig',1,'type','weibull');
    disp('Gender performance during CONTRAST task');
    doStaircase('threshold',stimulus.p.dualstaircase{2},'dispFig',1,'type','weibull');
    disp('Gender performance during SOLO task');
    doStaircase('threshold',usestair,'dispFig',1,'type','weibull');
catch
    warning('Peripheral task staircases failed');
end

%% Generate discrimination functions

% The idea here is just to compare the single and dual task performance on
% the same graphs.
gen_discFuncs(stimulus);

%% Generate Performance Plots

% The idea here is to have a plot that shows the dual task performance in
% comparison with the single task performance for both gender and
% contrast/noise at the same time.

gen_perf(stimulus);
%% To start off, let's get files

global analysis

[analysis.datFolder, analysis.anFolder] = getSubjDataFolder;

year = date;
year = year(end-1:end);

files = dir(sprintf('%s/%s*mat',analysis.datFolder,year));

%% Loop over files, fixing them for analysis

for fi = 1:length(files)
    curFile = sprintf('%s/%s',analysis.datFolder,files(fi).name);
    fixFile = sprintf('%s/%s',analysis.datFolder,['f' files(fi).name]);
    if isfile(fixFile)
        disp(sprintf('(noisecon_analysis) Found fixed file for %s, ignoring...',curFile));
    else
        success = genExp(curFile,fixFile);
        if success
            disp(sprintf('(noisecon_analysis) Output in: %s',fixFile));
        else
            disp(sprintf('(noisecon_analysis) File %s was not fixed!!!!',curFile));
        end
    end
end

%% Loop over files and load them into expHolder

expHolder = cell(length(files));
for fi = 1:length(files)
    fixFile = sprintf('%s/%s',analysis.datFolder,['f' files(fi).name]);
    load(fixFile);
    expHolder{fi} = exp;
end

%% ANALYSIS TIME

% LIST OF ANALYSES THAT I EXPECT TO BE IMPORTANT:
%   -   Do errors during a dual trial cause the trial to fail
%       a) Are errors during the first half vs second half more important?
%   -   Is reaction time different?
%       a) Male vs. Female vs. None peripheral
%   -   Are specific images problematic for the peripheral task

%% Check for image problems in peripheral task

imagePerf.M = cell(1,100); % not totally sure of actual number of images
imagePerf.F = cell(1,100);
type = {'M','F'};

tCounter = 0;

for i = 1:length(expHolder)
    cExp = expHolder{i};
    main = cExp{1};
    per = cExp{2};
    rVars = cExp{3}.runVars;
    
    for t = 1:per.nTrials
        if per.randVars.respond(t) == 1
            gen = type{per.randVars.tGen(t)};
            corr = per.correct(t);
            if isnan(corr), corr = 0; end
            imagePerf.(gen){per.randVars.tImage(t)}(end+1) = corr;
            tCounter = tCounter + 1;
        end
    end
end
for m = 1:length(imagePerf.M)
    imagePerf.Mu(m) = nanmean(imagePerf.M{m});
    imagePerf.MN(m) = length(imagePerf.M{m});
end
for f = 1:length(imagePerf.F)
    imagePerf.Fu(f) = nanmean(imagePerf.F{f});
    imagePerf.FN(f) = length(imagePerf.F{f});
end
% figure
% plot(imagePerf.Mu);

%% Find images that have zero performance && Save
zero_perf.M = find(logical(imagePerf.Mu<.1).*logical(imagePerf.MN>2));
zero_perf.F = find(logical(imagePerf.Fu<.1).*logical(imagePerf.FN>2));

imgFolderM = '~/proj/att_awe/images/brazil_faces/m/';
imgFolderF = '~/proj/att_awe/images/brazil_faces/f/';

filesM = dir(strcat(imgFolderM,'*jpg'));
filesF = dir(strcat(imgFolderF,'*jpg'));

filesF = filesF(zero_perf.F);
filesM = filesM(zero_perf.M);

zero_perf.filesF = filesF;
zero_perf.filesM = filesM;

zero_perf.trials = tCounter;

save(strcat(analysis.anFolder,'/zero_perf.mat'),'zero_perf');

%% Check for reaction time differences

%               type = M1 F2 O3
% using gender(single/dual,trial,order,correct,type) = RT
gender = gen_perRT(expHolder);

%% Display RT

disp('First Stim');
disp(sprintf('RT Gender: Single Task, Male %0.3f Female %0.3f None %0.3f',nanmean(nanmean(gender(1,:,1,:,1))), ...
    nanmean(nanmean(gender(1,:,1,:,2))), nanmean(nanmean(gender(1,:,1,:,3)))));
disp(sprintf('RT Gender: Dual Task, Male %0.3f Female %0.3f None %0.3f',nanmean(nanmean(gender(2,:,1,:,1))), ...
    nanmean(nanmean(gender(2,:,1,:,2))), nanmean(nanmean(gender(2,:,1,:,3)))));
disp('Second Stim');
disp(sprintf('RT Gender: Single Task, Male %0.3f Female %0.3f None %0.3f',nanmean(nanmean(gender(1,:,2,:,1))), ...
    nanmean(nanmean(gender(1,:,2,:,2))), nanmean(nanmean(gender(1,:,2,:,3)))));
disp(sprintf('RT Gender: Dual Task, Male %0.3f Female %0.3f None %0.3f',nanmean(nanmean(gender(2,:,2,:,1))), ...
    nanmean(nanmean(gender(2,:,2,:,2))), nanmean(nanmean(gender(2,:,2,:,3)))));

%% Check for errors during dual trials

for i = 5:length(expHolder)%1:length(expHolder)
    cExp = expHolder{i};
    main = cExp{1};
    per = cExp{2};
    rVars = cExp{3}.runVars;
    if rVars.dual == 1
        % This is a dual run that is loaded
        correct = main.randVars.interval==main.response;
        pcorrect = zeros(size(correct));
        for j = 1:length(correct)
            % find the corresponding trials in the peripheral task, check
            % whether these were correct
            pcorrect(j) = nansum(per.correct(per.randVars.mainTrialNum==j));
            RT{pcorrect(j)+1}(end+1) = main.reactionTime(j);
        end
        perf(1) = sum(correct(pcorrect==0))/sum(pcorrect==0);
        perf(2) = sum(correct(pcorrect==1))/sum(pcorrect==1);
        perf(3) = sum(correct(pcorrect==2))/sum(pcorrect==2);
    end
end

%% Print main RT influenced by dual

disp(sprintf('When no errors occurred, main task RT = %0.3f+-%0.3f, Perf = %0.2f%%',nanmean(RT{3}),1.96*nanstd(RT{3}),perf(3)));
disp(sprintf('When one error occurred, main task RT = %0.3f+-%0.3f, Perf = %0.2f%%',nanmean(RT{2}),1.96*nanstd(RT{2}),perf(2)));
disp(sprintf('When two errors occurred, main task RT = %0.3f+-%0.3f, Perf = %0.2f%%',nanmean(RT{1}),1.96*nanstd(RT{1}),perf(1)));

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

%% The following functions use the most recent stimulus file

load(fullfile(analysis.datFolder,files(end).name));

%% Generate discrimination functions

% The idea here is just to compare the single and dual task performance on
% the same graphs.
plotting = gen_discFuncs(stimulus);

%% Generate Performance Plots

% The idea here is to have a plot that shows the dual task performance in
% comparison with the single task performance for both gender and
% contrast/noise at the same time.

gen_perf(stimulus,plotting);
%% To start off, let's get files

if isempty(mglGetSID)
    error('Please set subject ID first.');
end

datFolder = sprintf('~/data/cue_noisecon/%s',mglGetSID);

year = date;
year = year(end-1:end);

files = dir(sprintf('%s/%s*mat',datFolder,year));

%% Loop over files, fixing them for analysis

exp_holder = {};
for fi = 1:length(files)
    curFile = sprintf('%s/%s',datFolder,files(fi).name);
    fixFile = sprintf('%s/%s',datFolder,['f' files(fi).name]);
    if isfile(fixFile)
        disp(sprintf('(noisecon_analysis) Found fixed file, loading... %s',fixFile));
        load(fixFile);
    else
        % No file, we need to fix the current one.
        disp(sprintf('(noisecon_analysis) Loading file %s',curFile));
        load(curFile);
        disp(sprintf('(noisecon_analysis) Fixing file %s',curFile));
        exp = getTaskParameters(myscreen,task);
        %% Peripheral Re-organization
        exp{2}.randVars.gender1 = zeros(size(exp{2}.randVars.position));
        exp{2}.randVars.gender2 = zeros(size(exp{2}.randVars.position));
        exp{2}.randVars.image1 = zeros(size(exp{2}.randVars.position));
        exp{2}.randVars.image2 = zeros(size(exp{2}.randVars.position));
        for i = 1:length(exp{2}.randVars.position)
            exp{2}.randVars.gender1(i) = exp{2}.randVars.gender{i}(1);
            exp{2}.randVars.gender2(i) = exp{2}.randVars.gender{i}(2);
            exp{2}.randVars.image1(i) = exp{2}.randVars.images{i}(2);
            exp{2}.randVars.image2(i) = exp{2}.randVars.images{i}(2);
        end
        exp{2}.randVars = rmfield(exp{2}.randVars,'gender');
        exp{2}.randVars = rmfield(exp{2}.randVars,'images');
        %% Main Re-organization
        exp{1}.randVars.image1 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.image2 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.image3 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.image4 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.gender1 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.gender2 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.gender3 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.gender4 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.contrast1_int1 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.contrast1_int2 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.contrast2_int1 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.contrast2_int2 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.contrast3_int1 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.contrast3_int2 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.contrast4_int1 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.contrast4_int2 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.noise1_int1 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.noise1_int2 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.noise2_int1 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.noise2_int2 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.noise3_int1 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.noise3_int2 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.noise4_int1 = zeros(size(exp{1}.randVars.deltaPed));
        exp{1}.randVars.noise4_int2 = zeros(size(exp{1}.randVars.deltaPed));
        for i = 1:length(exp{1}.randVars.deltaPed)
            exp{1}.randVars.image1(i) = exp{1}.randVars.imageNums{i}(1);
            exp{1}.randVars.image2(i) = exp{1}.randVars.imageNums{i}(2);
            exp{1}.randVars.image3(i) = exp{1}.randVars.imageNums{i}(3);
            exp{1}.randVars.image4(i) = exp{1}.randVars.imageNums{i}(4);
            exp{1}.randVars.gender1(i) = exp{1}.randVars.genderList{i}(1);
            exp{1}.randVars.gender2(i) = exp{1}.randVars.genderList{i}(2);
            exp{1}.randVars.gender3(i) = exp{1}.randVars.genderList{i}(3);
            exp{1}.randVars.gender4(i) = exp{1}.randVars.genderList{i}(4);
            exp{1}.randVars.contrast1_int1(i) = exp{1}.randVars.contrastList{i}(1,1);
            exp{1}.randVars.contrast1_int2(i) = exp{1}.randVars.contrastList{i}(1,2);
            exp{1}.randVars.contrast2_int1(i) = exp{1}.randVars.contrastList{i}(2,1);
            exp{1}.randVars.contrast2_int2(i) = exp{1}.randVars.contrastList{i}(2,2);
            exp{1}.randVars.contrast3_int1(i) = exp{1}.randVars.contrastList{i}(3,1);
            exp{1}.randVars.contrast3_int2(i) = exp{1}.randVars.contrastList{i}(3,2);
            exp{1}.randVars.contrast4_int1(i) = exp{1}.randVars.contrastList{i}(4,1);
            exp{1}.randVars.contrast4_int2(i) = exp{1}.randVars.contrastList{i}(4,2);
            exp{1}.randVars.noise1_int1(i) = exp{1}.randVars.noiseList{i}(1,1);
            exp{1}.randVars.noise1_int2(i) = exp{1}.randVars.noiseList{i}(1,2);
            exp{1}.randVars.noise2_int1(i) = exp{1}.randVars.noiseList{i}(2,1);
            exp{1}.randVars.noise2_int2(i) = exp{1}.randVars.noiseList{i}(2,2);
            exp{1}.randVars.noise3_int1(i) = exp{1}.randVars.noiseList{i}(3,1);
            exp{1}.randVars.noise3_int2(i) = exp{1}.randVars.noiseList{i}(3,2);
            exp{1}.randVars.noise4_int1(i) = exp{1}.randVars.noiseList{i}(4,1);
            exp{1}.randVars.noise4_int2(i) = exp{1}.randVars.noiseList{i}(4,2);
        end
        exp{1}.randVars = rmfield(exp{1}.randVars,'imageNums');
        exp{1}.randVars = rmfield(exp{1}.randVars,'genderList');
        exp{1}.randVars = rmfield(exp{1}.randVars,'contrastList');
        exp{1}.randVars = rmfield(exp{1}.randVars,'noiseList');
        %% Stimulus Fix:
        % Staircase is, (NOISE1/CON2),CUE1/CUE4,PED1:3
        exp{1}.stairNoise = squeeze(stimulus.staircase(1,:,:));
        exp{1}.stairNoiseDual = squeeze(stimulus.dualstaircase(1,:,:));
        exp{1}.stairCon = squeeze(stimulus.staircase(2,:,:));
        exp{1}.stairConDual = squeeze(stimulus.dualstaircase(2,:,:));
        exp{1}.p.stair = squeeze(stimulus.p.staircase);
        exp{1}.p.stairDual = squeeze(stimulus.p.dualstaircase);
        %% run Variables
        exp{3}.runVars.blocks = stimulus.blockList;
        exp{3}.runVars.dual = stimulus.dualList;
        exp{3}.runVars.pedestals = stimulus.pedestals;
        exp{3}.runVars.trainingRun = stimulus.training;
        %% Save
        disp(sprintf('(noisecon_analysis) Saving fix file %s',fixFile));
        exp_holder{stimulus.counter} = exp;
        save(fixFile,'exp');
    end
end

%% Generate staircase graphs
% this just generates the staircases for any existing runs

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

disp('Computing Weibull functions. CAUTION: Check all Weibull functions for accuracy');
check = 0;

dispTexts = {'Noise','Contrast'};
colorOpts = [1 0 0
             0 1 0];
typePs = {'noise','contrast'};
% plotting = zeros(2,3,2,2);
for num = 1:2
    dispText = dispTexts{num};
    color = colorOpts(num,:);
    typeP = typePs{num};
    
%     for cues = 1:2
%         for ped = 1:3
%             try
%                 if check
%                     out1 = doStaircase('threshold',stimulus.staircase{num,cues,ped},'dispFig',1,'type','weibull'); % noise, 1 cue, lowest
%                     out2 = doStaircase('threshold',stimulus.dualstaircase{num,cues,ped},'dispFig',1,'type','weibull'); % noise, 1 cue, lowest
%                     keyboard
%                 else
%                     out1 = doStaircase('threshold',stimulus.staircase{num,cues,ped},'type','weibull'); % noise, 1 cue, lowest
%                     out2 = doStaircase('threshold',stimulus.dualstaircase{num,cues,ped},'type','weibull'); % noise, 1 cue, lowest
%                 end
%                 plotting(cues,ped,1,num) = out1.threshold;
%                 plotting(cues,ped,2,num) = out2.threshold;
%             catch
%                 plotting(cues,ped,1,num) = -1;
%                 plotting(cues,ped,2,num) = -1;
%             end
%         end
%     end
    hold off
    % Discrimination function plots
    figure
    hold on
    title(sprintf('%s',dispText));
    plot(stimulus.pedestals.(typeP)(2:4),plotting(1,:,1,num),'-','Color',color);
    plot(stimulus.pedestals.(typeP)(2:4),plotting(2,:,1,num),'--','Color',color);
    plot(stimulus.pedestals.(typeP)(2:4),plotting(1,:,2,num),'-','Color',.5*color);
    plot(stimulus.pedestals.(typeP)(2:4),plotting(2,:,2,num),'--','Color',.5*color);
    if num == 1
        axis([stimulus.pedestals.(typeP)(2)-.1 stimulus.pedestals.(typeP)(4)+.1 0 .8]);
    else
        axis([stimulus.pedestals.(typeP)(2)-.05 stimulus.pedestals.(typeP)(4)+.05 0 .4]);
    end
    legend('Focal, Single Task','Distributed, Single Task','Focal, Dual Task','Distributed, Dual Task');
    print(gcf,'-dpdf',sprintf('~/proj/att_awe/analysis/figures/%sDiscriminationFunction',dispText));
    hold off
end

%% Generate Performance Plots

% The idea here is to have a plot that shows the dual task performance in
% comparison with the single task performance for both gender and
% contrast/noise at the same time.

% First let's choose what plot

for dual = 1:2
    if dual == 1
        stairtype = 'staircase';
        typeD = 'single';
    else
        stairtype = 'dualstaircase';
        typeD = 'dual';
    end
    for type = 1:2
        plotting = zeros(2,3);
        if type == 1
            % noise
            if isfield(stimulus.pedestals,'SnR')
                typeP = 'SnR';
            else
                typeP = 'noise';
            end
            dispText = 'Noise';
            num = 1;
        else
            typeP = 'contrast';
            dispText = 'Contrast';
            num = 2;
        end

        % Okay now we know what we're calculating for, so let's get the
        % performance

        % Get the main task performance
        for cues = 1:2
            for ped = 1:3
                try
                    out = doStaircase('threshold',stimulus.(stairtype){num,cues,ped},'type','weibull'); % noise, 1 cue, lowest
                    plotting(cues,ped) = out.threshold;
                catch
                    plotting(cues,ped) = -1;
                end
            end
        end

        main.(dispText).(typeD) = mean(plotting,2);
    end
end

% Main task performance
mainConFocalPerf = main.Contrast.single(1);
mainConDistPerf = main.Contrast.single(2);
mainNoiseFocalPerf = main.Noise.single(1);
mainNoiseDistPerf = main.Noise.single(2);

mainConFocalDualPerf = main.Contrast.dual(1);
mainConDistDualPerf = main.Contrast.dual(2);
mainNoiseFocalDualPerf = main.Noise.dual(1);
mainNoiseDistDualPerf = main.Noise.dual(2);

% Get the peripheral task performance
%%%% check for >2 task sets
if length(stimulus.p.dualstaircase{1}) > 1
%     genderNoisePerf = doStaircase('threshold',stimulus.p.dualstaircase{1}(3:end),'type','weibull');
    genderNoisePerf = doStaircase('threshold',stimulus.p.dualstaircase{1}(3:end));
else
%     genderNoisePerf = doStaircase('threshold',stimulus.p.dualstaircase{1},'type','weibull');
    genderNoisePerf = doStaircase('threshold',stimulus.p.dualstaircase{1});
end
if length(stimulus.p.dualstaircase{2}) > 1
    genderConPerf = doStaircase('threshold',stimulus.p.dualstaircase{2}(3:end),'type','weibull');
%     genderConPerf = doStaircase('threshold',stimulus.p.dualstaircase{2}(3:end));
else
%     genderConPerf = doStaircase('threshold',stimulus.p.dualstaircase{2},'type','weibull');
    genderConPerf = doStaircase('threshold',stimulus.p.dualstaircase{2});
end
% genderPerf = doStaircase('threshold',stimulus.p.staircase(2:end),'type','weibull');
genderPerf = doStaircase('threshold',stimulus.p.staircase(1:end-1));
% gNPerf = genderNoisePerf.threshold;
% gCPerf = genderConPerf.threshold;
% gPerf = genderPerf.threshold;
gNPerf = genderNoisePerf.meanOfLast7Reversals;
gCPerf = genderConPerf.meanOfLast7Reversals;
gPerf = genderPerf.meanOfLast7Reversals(end);
% % % % % % % % % Normalize
% % % % % % % % gNPerf_N = gNPerf;
% % % % % % % % gCPerf_N = gCPerf;
% % % % % % % % gPerf = gPerf;
% Plot
figure
hold on
title('Dual Task Performance');
% singles
plot(0,gPerf,'*r');
ylabel('Gender (SOA ms)');
xlabel('Contrast/Noise Performance (delta)');
plot(mainNoiseFocalPerf,0,'*g');
text(mainNoiseFocalPerf,.01,sprintf('%0.2f',mainNoiseFocalPerf/(1-mainNoiseFocalPerf)));
plot(mainNoiseDistPerf,0,'*c');
text(mainNoiseDistPerf,.01,sprintf('%0.2f',mainNoiseDistPerf/(1-mainNoiseDistPerf)));
plot(mainConFocalPerf,0,'*r');
plot(mainConDistPerf,0,'*m');
% duals
plot(mainNoiseFocalDualPerf,gNPerf,'*g');
text(mainNoiseFocalDualPerf,gNPerf+.01,sprintf('%0.2f',mainNoiseFocalDualPerf/(1-mainNoiseFocalDualPerf)));
plot(mainNoiseDistDualPerf,gNPerf,'*c');
text(mainNoiseDistDualPerf,gNPerf+.01,sprintf('%0.2f',mainNoiseDistDualPerf/(1-mainNoiseDistDualPerf)));
plot(mainConFocalDualPerf,gCPerf,'*r');
plot(mainConDistDualPerf,gCPerf,'*m');
% lines
plot(0:mainNoiseDistPerf/10:mainNoiseDistPerf,repmat(gPerf,1,11),'--r');
plot(0:.25:.5,repmat(.25,1,3),'--k');
plot(repmat(mainNoiseFocalPerf,1,11),0:gPerf/10:gPerf,'--g');
plot(repmat(mainNoiseDistPerf,1,11),0:gPerf/10:gPerf,'--c');
plot(repmat(mainConFocalPerf,1,11),0:gPerf/10:gPerf,'--r');
plot(repmat(mainConDistPerf,1,11),0:gPerf/10:gPerf,'--m');
plot(repmat(.5,1,3),0:.125:.25,'--k');
% legend
% legend('Gender','Focal Noise','Dist Noise','Focal Contrast','Dist Contrast');
% distance lines
xp = [mainConFocalPerf mainConFocalDualPerf];
yp = [gPerf gCPerf];
distCF = sqrt((xp(2)-xp(1))^2 + (yp(2)-yp(1))^2);
plot(xp,yp,'-k');
text(mean(xp),mean(yp),sprintf('Dist: %1.2f',distCF));

xp = [mainNoiseFocalPerf mainNoiseFocalDualPerf];
yp = [gPerf gNPerf];
distNF = sqrt((xp(2)-xp(1))^2 + (yp(2)-yp(1))^2);
plot(xp,yp,'-k');
text(mean(xp),mean(yp),sprintf('Dist: %1.2f',distNF));

xp = [mainConDistPerf mainConDistDualPerf];
yp = [gPerf gCPerf];
distCD = sqrt((xp(2)-xp(1))^2 + (yp(2)-yp(1))^2);
plot(xp,yp,'-k');
text(mean(xp),mean(yp),sprintf('Dist: %1.2f',distCD));

xp = [mainNoiseDistPerf mainNoiseDistDualPerf];
yp = [gPerf gNPerf];
distND = sqrt((xp(2)-xp(1))^2 + (yp(2)-yp(1))^2);
plot(xp,yp,'-k');
text(mean(xp),mean(yp),sprintf('Dist: %1.2f',distND));

disp(sprintf('Euclidian Distance for Contrast: %1.2f, for Noise: %1.2f',sqrt(distCF^2+distCD^2),sqrt(distNF^2+distND^2)));
disp('This distance measurement is innacurate :), noise should be on a log scale and euclidian distance isn''t really correct here');
 hold off
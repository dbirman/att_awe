% cue_fade
%
%      usage: myscreen=cue_fade()
%         by: daniel birman
%       date: 11/10/14
%    purpose: contrast change detection with cued selective attention.
%
%        use: call cue_fade() to initialize. The first
%             run takes significantly longer due to loading stimuli.
%
%      flags: stimFileNum (-1/#) - Load a specific stimfile from a
%             subject's folder. Defaults to the last stimfile. Warning:
%             Only the first stimfile is saved with the image file data,
%             subsequent stimfiles only contain run data.
%             plots (0/1) - Displays staircase plots (and estimated
%             psychophysic functions)
%             practice (def=0) - Displays on ten trials
%

function [myscreen] = cue_fade(varargin)

global stimulus
%% Initialize Variables

% add arguments later
widthDeg = []; heightDeg = [];
stimFileNum = [];
plots = [];
practice = [];
getArgs(varargin,{'widthDeg=5.5', 'heightDeg=5.5', ...
    'stimFileNum=-1', ...
    'dual=0','plots=1','practice=0'});

stimulus.counter = 1; % This keeps track of what "run" we are on.
%% Setup Screen

screen.screenNumber = 2;
myscreen = initScreen(screen);

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/cue_fade/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/cue_fade/%s/1*mat',mglGetSID));

    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/cue_fade/%s/%s',mglGetSID,fname));
        stimulus.staircase = s.stimulus.staircase;
        stimulus.counter = s.stimulus.counter + 1;

        % load blocks too
        stimulus.blocks = s.stimulus.blocks;
        stimulus.blocks.loaded = 1;

        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(fade) Data file: %s loaded, this is run #%i',fname,stimulus.counter));
    end
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

stimulus.responseKeys = [10 9]; % corresponds to CHANGE - NO CHANGE

stimulus.colors.rmed = 127.5;
stimulus.colors.mrmax = 255;
stimulus.colors.mrmin = 0;

stimulus.basepdf = normpdf(stimulus.colors.mrmin:stimulus.colors.mrmax,stimulus.colors.rmed,75);


%% MGL Parameters
mglTextSet('Helvetica',32,255,... % this doesn't work... not sure why
0,0,0,0,0,0,0);
    
%% Initialize Images

% set initial thresholds
stimulus.pedestals.contrast = [ .4 .6 .9];
stimulus.pedestals.noise = .4;
stimulus.baseThresh(1) = .35;

stimulus.nPedestals = length(stimulus.pedestals.contrast);

% load images
stimulus.widthDeg = widthDeg;
stimulus.heightDeg = heightDeg;
stimulus.posx = [-3.8 +3.8 +3.8 -3.8];
stimulus.posy = [+3.8 +3.8 -3.8 -3.8];
categories = {'m' 'f'};
name = getenv('USER');
stimulus.widthPix = 400;
stimulus.heightPix = 400;
stimulus.imageDirMain = fullfile(sprintf('/Users/%s/proj/att_awe/images/all_faces2/',name));
dispLoadFig = 0; keepAspectRatio = 0;

stimulus = InitStimFade(stimulus,categories,dispLoadFig,keepAspectRatio);

%% Initialize Mask

gaussianWin = mglMakeSmoothBorderMask(1,1,.9,.1,.1,399,399);
win = 255-255*gaussianWin;
% win = 255-255*(gaussianWin>0);
stimulus.mask = ones(size(win,1),size(win,2),3)*stimulus.colors.rmed;
stimulus.mask(:,:,4) = win;
stimulus.maskTex = mglCreateTexture(stimulus.mask);

%% Choose step sizes for changes

stimulus.steps = 20;

%% Setup Task

% This is the contrast change detection task
task{1}{1}.waitForBacktick = 1;
stimulus.seg.ITI = 1; % the ITI is either 20s (first time) or 1s
stimulus.seg.cue = 2; % the cue is on for 1s
stimulus.seg.stim_1hold = 3; % the stimulus is on for 1s
stimulus.seg.stim_2chng = 4;
stimulus.seg.ISI = 5;
stimulus.seg.resp = 6;
task{1}{1}.segmin = [1 .8 .1 .5 .25 1.3];
task{1}{1}.segmax = [1 .8 .5 .5 .75 1.3];
task{1}{1}.synchToVol = [0 0 0 0 0 0];
task{1}{1}.getResponse = [0 0 0 0 0 1];
task{1}{1}.parameter.blockTrialNum = 1:20; % we just need this to have the right number of trials in each block, we will add our own parameters at each trialstart
task{1}{1}.numBlocks = 5;

%% Block variables

% Within a "block" we are going to show all of the possible target (1:4) x
% attention conditions (cue 1/4) twice. This is a total of sixteen trials.
% These values are randomized for each block and are pulled as a random
% permutation of a matrix stored in stimulus.blocks.curTrialPerm

task{1}{1}.randVars.calculated.target = nan; % This is the target that gets cued
task{1}{1}.randVars.calculated.changeTarget = nan; % This is the actual location that changes
task{1}{1}.randVars.calculated.cues = nan;
task{1}{1}.randVars.calculated.change = nan; % will the stimulus actually change
task{1}{1}.randVars.calculated.catch = 0;
task{1}{1}.randVars.calculated.catchTask = nan;

% Each block has a very specific set of stimuli which are maintained across
% the block. Again, these are counterbalanced across blocks in a random
% permutation. There are 3 contrast pedestals, 4 locations, 2 image genders

%%
if practice
    task{1}{1}.numBlocks = 1;
end

%% Tracking

% these are variables that we want to track for later analysis.

task{1}{1}.randVars.calculated.pedestalList = nan(4);
task{1}{1}.randVars.calculated.genderList = nan(4);
task{1}{1}.randVars.calculated.maxContrast = nan;
task{1}{1}.randVars.calculated.deltaPed = nan;
task{1}{1}.randVars.calculated.imageNums = nan(1,4);
task{1}{1}.randVars.calculated.correct = nan;
task{1}{1}.randVars.calculated.trialNum = nan;

stimulus.curTrial = 0;

%% Block setup
if isfield(stimulus,'blocks') && isfield(stimulus.blocks,'loaded')
    % We already have our blocks
    stimulus.blocks = rmfield(stimulus.blocks,'loaded'); % remove the load field, otherwise it gets saved across runs
    % Make sure we have enough pedestalListOptions
    if stimulus.blocks.counter + 7 > size(stimulus.blocks.pedestalListOptions,1)
        % Not enough, let's double the size of pedestalListOptions
        stimulus.blocks.pedestalListOptions = [stimulus.blocks.pedestalListOptions ; ...
            stimulus.blocks.pedestalListOptions(randperm(size(stimulus.blocks.pedestalListOptions,1)),:)];
    end
else
    % This is the first run, build up the blocks.
    
    % In the permutation matrix we want to track things that are important
    % within a single block and need to be counterbalanced. Target #, cues,
    % and change. [Targets;Cues;Change]
    stimulus.blocks.permutationMatrix = [repmat(1:4,1,4)',repmat([1,1,1,1,2,2,2,2],1,2)',[zeros(1,8),ones(1,8)]'];
    stimulus.blocks.permutationMatrix = [stimulus.blocks.permutationMatrix ; 1 1 1; 2 1 1; 3 1 1; 4 1 1];
    stimulus.blocks.permutationMatrix = [stimulus.blocks.permutationMatrix ; 1 1 0; 2 1 0; 3 1 0; 4 1 0];
    % Now we want to set up the target pedestals we'll use for each of the
    % blocks (stored in thistrial.pedestalList). On each block we'll use on
    % of these permutations, until we run out and we reset the list.
    counter = 1;
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    stimulus.blocks.pedestalListOptions(counter,:) = [i,k,k,l];
                    counter = counter + 1;
                end
            end
        end
    end
    stimulus.blocks.pedestalListOptions = stimulus.blocks.pedestalListOptions(randperm(size(stimulus.blocks.pedestalListOptions,1)),:);
    stimulus.blocks.counter = 0;
end

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],@startBlockCallback);
end

mglClearScreen(stimulus.colors.rmed/255);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimulus.initStair
    % We are starting our staircases from scratch
    disp(sprintf('(fade) Initializing staircases'));
    stimulus = initStaircase(stimulus);
else
    disp('(fade) Re-using staircase from previous run...');
    % Reset staircase if necessary
    checkStaircaseStop(stimulus);
end

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%% Main Task Loop

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

% task ended
mglClearScreen(stimulus.colors.rmed/255);
mglTextDraw('Run complete... please wait.',[0 0]);
mglFlush
myscreen.flushMode = 1;
disp('(fade) Run ending...');

if plots
    disp('(fade) Displaying plots');
    dispStaircase(stimulus);
end

% delete texture
if isfield(stimulus,'flyTex')
    for x = 1:size(stimulus.flyTex,1)
        mglDeleteTexture(stimulus.flyTex{x});
    end
end
mglDeleteTexture(stimulus.maskTex);
stimulus = rmfield(stimulus,'flyTex');

%% Save file (not, first file will always have full textures)
if stimulus.counter > 1
    % Temporarily remove raw files to limit stimfile size, then replace these
    % so they can be used on the next run.
    stimbackup = stimulus;
    stimulus = rmfield(stimulus,'raw');
    stimulus = rmfield(stimulus,'images');
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

if stimulus.counter > 1
    stimulus = stimbackup;
    clear stimbackup
end

%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT OVER: HELPER FUNCTIONS FOLLOW %%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Trial %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startTrialCallback(task,myscreen)

myscreen.flushMode = 0;
% We need to do a number of things here. First we need to set up the
% missing thistrial parameters. Second, we will pre-calculate the textures
% that will be displayed on this particular trial. 

global stimulus

stimulus.curTrial = stimulus.curTrial + 1;
task.thistrial.trialNum = stimulus.curTrial;

task.thistrial.target = stimulus.blocks.curTrialPerm(task.thistrial.blockTrialNum,1);
task.thistrial.cues = stimulus.blocks.curTrialPerm(task.thistrial.blockTrialNum,2);
task.thistrial.change = stimulus.blocks.curTrialPerm(task.thistrial.blockTrialNum,3);
if task.thistrial.cues == 1
    % only on focal trials
    if rand < .2
        locs = [1 2 3 4];
        locs = locs(locs~=task.thistrial.target);
        task.thistrial.changeTarget = locs(randi(3));
        task.thistrial.catch = 1;
        task.thistrial.catchTask = randi(2);
        % If this is a catch trial, extend the response window to 2.5 s
        task.thistrial.seglen(stimulus.seg.resp) = 2.5;
    else
        task.thistrial.changeTarget = task.thistrial.target;
    end
else
    task.thistrial.changeTarget = task.thistrial.target;
end
task.thistrial.pedestalList = stimulus.blocks.curpedestalList;
task.thistrial.genderList = stimulus.blocks.curGenderList;
task.thistrial.maxContrast = stimulus.blocks.curMaxContrast;

% Get Delta
[task.thistrial.deltaPed, stimulus] = getDeltaPed(task,stimulus, ...
    task.thistrial.pedestalList(task.thistrial.target),task.thistrial.cues);

catchType = {'regular','catch'};
changeType = {'no change','change'};
% Display info
trialType = {'contrast','noise'};
if isnan(task.thistrial.catchTask)
    type = 1;
else
    type = task.thistrial.catchTask;
end
disp(sprintf('(fade) Trial %i is a %s trial. Displaying with %0.2f %s - %0.2f delta. With %s',task.thistrial.trialNum,catchType{task.thistrial.catch+1}, ...
    stimulus.pedestals.(trialType{type})(task.thistrial.pedestalList(task.thistrial.changeTarget)),trialType{type},task.thistrial.deltaPed,...
    changeType{task.thistrial.change+1}));

% Build changeTex
if task.thistrial.change
    stimulus = buildChangeTex(task,stimulus);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus

switch task.thistrial.thisseg
    case stimulus.seg.ITI
        stimulus.live.fixColor = 0;
        stimulus.live.changing = 0;
        stimulus.live.cues = 0;
        stimulus.live.faces = 0;
    case stimulus.seg.cue
        stimulus.live.fixColor = 0;
        stimulus.live.changing = 0;
        stimulus.live.cues = 1;
        stimulus.live.faces = 0;
    case stimulus.seg.stim_1hold
        stimulus.live.fixColor = 0;
        stimulus.live.changing = 0;
        stimulus.live.cues = 1;
        stimulus.live.faces = 1;
    case stimulus.seg.stim_2chng
        stimulus.live.fixColor = 0;
        stimulus.live.changing = task.thistrial.change;
        stimulus.live.cues = 1;
        stimulus.live.faces = 1;
    case stimulus.seg.ISI
        stimulus.live.fixColor = 0;
        stimulus.live.changing = 0;
        stimulus.live.cues = 0;
        stimulus.live.faces = 0;
    case stimulus.seg.resp
        stimulus.live.fixColor = 1;
        stimulus.live.changing = 0;
        stimulus.live.cues = 1;
        stimulus.live.faces = 0;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
global stimulus

mglClearScreen(stimulus.colors.rmed/255);

upFix(stimulus);
if stimulus.live.cues
    upCues(task,stimulus);
end
if stimulus.live.faces
    upFaces(task,stimulus);
end

%%
function upFix(stimulus)

mglFixationCross(1,1,stimulus.live.fixColor);

%%
function upCues(task,stimulus)

if task.thistrial.thisseg == stimulus.seg.resp
    tLoc = task.thistrial.changeTarget;
else
    tLoc = task.thistrial.target;
end

ang = [180 0 0 180];
if task.thistrial.cues == 1 || task.thistrial.thisseg == stimulus.seg.resp
    % We can't just display all of the lines, we just want one line
    usePos1 = stimulus.posx(tLoc);
    usePos2 = stimulus.posy(tLoc);
    ang = ang(tLoc);
else
    usePos1 = stimulus.posx;
    usePos2 = stimulus.posy;
end
ang = d2r(ang);
mglLines2(cos(ang + atan(usePos1./usePos2))*.1,sin(ang + atan(usePos1./usePos2))*.1, ...
    cos(ang + atan(usePos1./usePos2))*.75,sin(ang + atan(usePos1./usePos2))*.75,1,stimulus.live.fixColor);

function upFaces(task,stimulus)


for imagePos = 1:4
    if stimulus.live.changing && imagePos == task.thistrial.changeTarget
        time = task.thistrial.seglen(stimulus.seg.stim_2chng)*1000;
        rNum = length(stimulus.changeTex);
        % figure out how long it's been since this segment started
        changeTime = (mglGetSecs - task.thistrial.segStartSeconds) * 1000;
        % down only code
        frame = round(changeTime/(time/(rNum-1)))+1;
        if frame < 1, frame = 1; end
        if frame > rNum, frame = rNum; end
        curimage = stimulus.changeTex{frame};
    else
        % Get the image from the textures built on the fly
        curimage = stimulus.flyTex{imagePos};
    end
    % Push image to buffer
    mglBltTexture(curimage,[stimulus.posx(imagePos) stimulus.posy(imagePos) stimulus.widthDeg stimulus.heightDeg]);
    mglBltTexture(stimulus.maskTex,[stimulus.posx(imagePos) stimulus.posy(imagePos) stimulus.widthDeg stimulus.heightDeg]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus

responseText = {'Incorrect','Correct'};
fixColors = {[1 0 0],[0 1 0]};

if any(task.thistrial.whichButton == stimulus.responseKeys)
    if task.thistrial.gotResponse == 0
        % Store whether this was correct
        task.thistrial.correct = find(stimulus.responseKeys==task.thistrial.whichButton) == task.thistrial.change+1;
        stimulus.live.fixColor = fixColors{task.thistrial.correct+1};
        disp(sprintf('(fade) Response %s',responseText{task.thistrial.correct+1}));
        if ~task.thistrial.catch
            stimulus.staircase{task.thistrial.cues,task.thistrial.pedestalList(task.thistrial.target)} = ...
                doStaircase('update',stimulus.staircase{task.thistrial.cues,task.thistrial.pedestalList(task.thistrial.target)},task.thistrial.correct);
        else
            stimulus.stairCatch = doStaircase('update',stimulus.stairCatch,task.thistrial.correct);
        end
    else
        disp(sprintf('(fade) Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Block Call Back %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = startBlockCallback(task, myscreen)

global stimulus


% clear screen
mglClearScreen(stimulus.colors.rmed/255);
if task.trialnum==1
    mglTextDraw('Run Starting...',[0 0]);
end
mglFlush

myscreen.flushMode = 1;
% Increment block (it starts at 0)
stimulus.blocks.counter = stimulus.blocks.counter + 1;


% let the user know
disp(sprintf('(fade) Starting block number: %i',stimulus.blocks.counter));

stimulus.blocks.curTrialPerm = stimulus.blocks.permutationMatrix(randperm(size(stimulus.blocks.permutationMatrix,1)),:);
stimulus.blocks.curpedestalList = stimulus.blocks.pedestalListOptions(stimulus.blocks.counter,:);
stimulus.blocks.curMaxContrast = 1;
% stimulus.blocks.curMaxContrast = max(stimulus.pedestals.contrast(stimulus.blocks.curpedestalList));

% setGammaTableForMaxContrast(stimulus.blocks.curMaxContrast);

% now that we know the current contrasts values we can build the flyTex
% values for this current trial
stimulus.blocks.curGenderList = randi(2,4);
stimulus.blocks.curImageList = [-1 -1 -1 -1];
for i = 1:4
    stimulus.blocks.curPMask(i,:) = rand(1,length(stimulus.raw{1}.halfFourier{1}.phase))*2*pi;
end
while length(unique(stimulus.blocks.curImageList))~=4
    for img = 1:4
        iNum = randi(size(stimulus.images{stimulus.blocks.curGenderList(img)},3));
        used(img) = iNum;
        stimulus.blocks.curImageList(img) = iNum;
    end
end

stimulus = getBlockImages(stimulus);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              HELPER FUNCTIONS                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stimulus = buildChangeTex(task,stimulus)

% We want to build a distribution of textures that will be used across 200
% ms. We will make the assumption that the current screen refresh rate is
% 150 hz, so we need 15 images down. Note since this is trial start (the 
% ITI) it doesn't matter if we drop frames right now.

highValue = stimulus.pedestals.contrast(task.thistrial.pedestalList(task.thistrial.changeTarget));

lowValue = highValue - task.thistrial.deltaPed;

if lowValue <= 0
    disp('(fade) Required contrast range dropped below zero, truncating');
    lowValue = .01;
end

if task.thistrial.catchTask == 1
    range = fliplr(lowValue : (highValue-lowValue)/(stimulus.steps-1) : highValue);
else
    range = fliplr(0:task.thistrial.deltaPed/(stimulus.steps-1):task.thistrial.deltaPed);
end

if isfield(stimulus,'changeTex')
    for i = 1:length(stimulus.changeTex)
        mglDeleteTexture(stimulus.changeTex{i});
    end
else
    stimulus.changeTex = cell(1,length(range));
end

    % get the image
img = stimulus.images{stimulus.blocks.curGenderList(task.thistrial.changeTarget)}(:,:,stimulus.blocks.curImageList(task.thistrial.changeTarget));
nimg_f = getHalfFourier(img);
nimg_f.phase = stimulus.blocks.curPMask(task.thistrial.changeTarget,:);
noise = reconstructFromHalfFourier(nimg_f);
L0 = mean2(img);
if task.thistrial.catchTask == 2
    curCon = highValue;
    npdf = scalePdf(stimulus.basepdf,stimulus.colors.mrmin:stimulus.colors.mrmax,curCon / stimulus.blocks.curMaxContrast);
end
for i = 1:length(range)

%     else
tic
    if task.thistrial.catchTask == 2
        K = .4+range(i);
    else
        curCon = range(i);
        npdf = scalePdf(stimulus.basepdf,stimulus.colors.mrmin:stimulus.colors.mrmax,curCon / stimulus.blocks.curMaxContrast);
        K = .4;
    end

    image = L0 + sqrt(K) * (img - L0) + sqrt(1-K) * (noise - L0);
    % set the contrast
    % scale by the ratio
    % change the image to match the PDF
    imagef = (stimulus.colors.mrmax-stimulus.colors.mrmin)*histeq(image/255,npdf) + stimulus.colors.mrmin;
    disp(sprintf('time: %.3f',toc))
    tic
    stimulus.changeTex{i} = mglCreateTexture(imagef);
    disp(toc)
end

%% getDeltaPed

function [deltaPed, stimulus] = getDeltaPed(task,stimulus,tCon,cues)
if task.thistrial.catch
    [deltaPed, stimulus.staircatch] = doStaircase('testValue',stimulus.stairCatch);
else
    [deltaPed, stimulus.staircase{cues,tCon}] = doStaircase('testValue',stimulus.staircase{cues,tCon});
end


%%%%%%%%%%%%%%%%%%%%%%%%
%    getBlockImages    %
%%%%%%%%%%%%%%%%%%%%%%%%

function stimulus = getBlockImages(stimulus)

if isfield(stimulus,'flyTex')
    for i = 1:4
        mglDeleteTexture(stimulus.flyTex{i});
    end
end
for i = 1:4
    
    % get the image
    img = stimulus.images{stimulus.blocks.curGenderList(i)}(:,:,stimulus.blocks.curImageList(i));
    % Contrast
    
    curNoi = stimulus.pedestals.noise;
    
    img = pinkNoise(img,[],curNoi,stimulus.blocks.curPMask(i,:));
    
    % set the contrast
    curCon = stimulus.pedestals.contrast(stimulus.blocks.curpedestalList(i));
    % scale by the ratio
    npdf = scalePdf(stimulus.basepdf,stimulus.colors.mrmin:stimulus.colors.mrmax,curCon / stimulus.blocks.curMaxContrast);
    % change the image to match the PDF
    img = (stimulus.colors.mrmax-stimulus.colors.mrmin)*histeq(img/255,npdf) + stimulus.colors.mrmin;
    
    % add to flyTex
    stimulus.flyTex{i} = mglCreateTexture(img);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)

stimulus.stairCatch = doStaircase('init','fixed',...
    'fixedVals',[.075 .125 .175 .25]);
stimulus.staircase = cell(2,length(stimulus.pedestals.contrast));
stimulus.staircase{1,1} = doStaircase('init','fixed',...
    'fixedVals',[.05 .075 .09 .11 .13 .15 .175]);
stimulus.staircase{2,1} = doStaircase('init','fixed',...
    'fixedVals',[.075 .1 .12 .14 .16 .18 .205]);
for i = 2:length(stimulus.pedestals.contrast)
    stimulus.staircase{1,i} = stimulus.staircase{1,1};
    stimulus.staircase{2,i} = stimulus.staircase{2,1};
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus)

try
    taskOpts = {'contrast'};
    
    figure % this is the 'staircase' figure
    title('%s, Staircase plot (R->G->B high)');
    hold on
    drawing = {'-r' '-g' '-b'
                '--r' '--g' '--b'};
    for cues = 1:2
        for ped = 1:3
            try
                testV = [];
                for i = 1:length(stimulus.staircase{cues,ped})
                    testV = [testV stimulus.staircase{cues,ped}(i).testValues];
                end
                plot(testV,drawing{cues,ped});
            catch
            end
            try
                out = doStaircase('threshold',stimulus.staircase{cues,ped},'type','weibull'); % noise, 1 cue, lowest
                plotting(cues,ped) = out.threshold;
            catch
                plotting(cues,ped) = -1;
            end
        end
    end
    hold off
    figure
    hold on
    title(sprintf('%s, R->G->B High'));
    plot(stimulus.pedestals.(taskOpts{1})(1:3),plotting(1,:),'-r');
    plot(stimulus.pedestals.(taskOpts{1})(1:3),plotting(2,:),'--r');
    axis([stimulus.pedestals.(taskOpts{1})(1) stimulus.pedestals.(taskOpts{1})(3) 0 1]);
    hold off

catch
    disp('(fade) Figures were not generated successfully.');
end

%% checkStaircaseStop
function checkStaircaseStop(stimulus)
% Check both staircases
for cues = 1:2
    for ped = 1:3
        s = stimulus.staircase{cues,ped,stimulus.blocks.curTask};
        if doStaircase('stop',s)
            % this is a bit of a pain... you can't pass an initialThreshold
            % argument do doStaircase('init',s, ...), it ignores everything and
            % resets using the calculated threshold. Because you can't override it
            [args, vals, ~] = getArgs(s(1).initArgs);
            threshPos = -1;
            stepPos = -1;
            for i = 1:length(args)
                switch args{i}
                    case 'initialThreshold'
                        threshPos = i;
                    case 'initialStepsize'
                        stepPos = i;
                end
            end
            out = doStaircase('threshold',s);
            in = input(sprintf('Resetting Staircase... Estimate is: %1.2f. Reset ([Y]/[C]ustom/[O]riginal): ',out.threshold),'s');
            switch in
                case 'Y'
                    vals{threshPos} = out.threshold;
                    vals{stepPos} = out.threshold / 3;
                case 'C'
                    disp('Original values:');
                    disp(sprintf('%s: %0.2f',args{threshPos},num2str(vals{threshPos})));
                    val = str2double(input('New threshold value: ','s'));
                    vals{threshPos} = val;
                    vals{stepPos} = val / 3;
                case 'O'
            end
            disp('THIS CODE IS NOT CERTAIN TO WORK! CHECK THE OUTPUT!');
            stimulus.staircase{cues,ped,stimulus.blocks.curTask} = doStaircase('init',s,'initialThreshold',vals{threshPos},'initialStepsize',vals{stepPos});
        end
    end
end

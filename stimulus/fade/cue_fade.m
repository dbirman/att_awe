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
    
%% Gamma Table Initialization

% set the reserved colors
% stimulus.gammaTable(1:size(stimulus.colors.reservedColors,1),1:size(stimulus.colors.reservedColors,2))=stimulus.colors.reservedColors;
% 
% % get gamma table
% if ~isfield(myscreen,'gammaTable')
%   stimulus.linearizedGammaTable = mglGetGammaTable;
%   disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
%   disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
%   disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
%   disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
% end
% stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

%% Initialize Images

% set initial thresholds
stimulus.pedestals.contrast = [ .4 .6 .9];
stimulus.baseThresh(1) = .35;
% These noise levels correspond to an SnR of 
% noisevals = [1.75 1.25 .75 .25 -.25];
% stimulus.pedestals.noise = 1./(1+exp(noisevals));
% stimulus.pedestals.SnR = stimulus.pedestals.noise ./ (1-stimulus.pedestals.noise);
% baseThresh(1) = .55;
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
stimulus.seg.resp = 5;
task{1}{1}.segmin = [1 1 .1 .5 1.2];
task{1}{1}.segmax = [1 1 .5 .5 1.2];
task{1}{1}.synchToVol = [0 0 0 0 0];
task{1}{1}.getResponse = [0 0 0 1 1];
task{1}{1}.parameter.blockTrialNum = 1:16; % we just need this to have the right number of trials in each block, we will add our own parameters at each trialstart
task{1}{1}.numBlocks = 7;

%% Block variables

% Within a "block" we are going to show all of the possible target (1:4) x
% attention conditions (cue 1/4) twice. This is a total of sixteen trials.
% These values are randomized for each block and are pulled as a random
% permutation of a matrix stored in stimulus.blocks.curTrialPerm

task{1}{1}.randVars.calculated.target = nan; % This is the target that gets cued
task{1}{1}.randVars.calculated.changeTarget = nan; % This is the actual location that changes
task{1}{1}.randVars.calculated.cues = nan;
task{1}{1}.randVars.calculated.change = nan; % will the stimulus actually change

% Each block has a very specific set of stimuli which are maintained across
% the block. Again, these are counterbalanced across blocks in a random
% permutation. There are 3 contrast pedestals, 4 locations, 2 image genders

%%
if practice
    task{1}{1}.numBlocks = 1;
end

%% Tracking

% these are variables that we want to track for later analysis.

% task{1}{1}.randVars.calculated.noiseList = nan(4,2);
task{1}{1}.randVars.calculated.contrastList = nan(4);
task{1}{1}.randVars.calculated.maxContrast = nan;
task{1}{1}.randVars.calculated.deltaPed = nan;
task{1}{1}.randVars.calculated.imageNums = nan(1,4);
task{1}{1}.randVars.calculated.correct = nan;

%% Block setup
if isfield(stimulus,'blocks') && isfield(stimulus.blocks,'loaded')
    % We already have our blocks
    stimulus.blocks = rmfield(stimulus.blocks,'loaded'); % remove the load field, otherwise it gets saved across runs
    % Make sure we have enough contrastListOptions
    if stimulus.blocks.counter + 7 > size(stimulus.blocks.contrastListOptions,1)
        % Not enough, let's double the size of contrastListOptions
        stimulus.blocks.contrastListOptions = [stimulus.blocks.contrastListOptions ; ...
            stimulus.blocks.contrastListOptions(randperm(size(stimulus.blocks.contrastListOptions,1)),:)];
    end
else
    % This is the first run, build up the blocks.
    
    % In the permutation matrix we want to track things that are important
    % within a single block and need to be counterbalanced. Target #, cues,
    % and change. [Targets;Cues;Change]
    stimulus.blocks.permutationMatrix = [repmat(1:4,1,4)',repmat([1,1,1,1,2,2,2,2],1,2)',[zeros(1,8),ones(1,8)]'];
    % Now we want to set up the target pedestals we'll use for each of the
    % blocks (stored in thistrial.contrastList). On each block we'll use on
    % of these permutations, until we run out and we reset the list.
    counter = 1;
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    stimulus.blocks.contrastListOptions(counter,:) = [i,k,k,l];
                    counter = counter + 1;
                end
            end
        end
    end
    stimulus.blocks.contrastListOptions = stimulus.blocks.contrastListOptions(randperm(size(stimulus.blocks.contrastListOptions,1)),:);
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

% This code can be used to adjust the length of the ITI (or any other
% segment)
% if task.thistrial.blockTrialNum == 1
%     task.thistrial.seglen(1) = 1.5;
% else
%     task.thistrial.seglen(1) = 1;
% end

task.thistrial.target = stimulus.blocks.curTrialPerm(task.thistrial.blockTrialNum,1);
task.thistrial.cues = stimulus.blocks.curTrialPerm(task.thistrial.blockTrialNum,2);
task.thistrial.change = stimulus.blocks.curTrialPerm(task.thistrial.blockTrialNum,3);
if task.thistrial.cues == 1
    % only on focal trials
    if rand < .1
        locs = [1 2 3 4];
        locs = locs(locs~=task.thistrial.target);
        task.thistrial.changeTarget = locs(randi(3));
    else
        task.thistrial.changeTarget = task.thistrial.target;
    end
end
task.thistrial.contrastList = stimulus.blocks.curContrastList;
task.thistrial.maxContrast = stimulus.blocks.curMaxContrast;

disp(sprintf('(cue_fade) Starting trial %i',task.thistrial.blockTrialNum));

% Get Delta
[task.thistrial.deltaPed, stimulus] = getDeltaPed(stimulus, ...
    task.thistrial.contrastList(task.thistrial.target),task.thistrial.cues);

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
        stimulus.live.faces = 0;
    case stimulus.seg.cue
        stimulus.live.fixColor = 0;
        stimulus.live.changing = 0;
        stimulus.live.faces = 0;
    case stimulus.seg.stim_1hold
        stimulus.live.fixColor = 0;
        stimulus.live.changing = 0;
        stimulus.live.faces = 1;
    case stimulus.seg.stim_2chng
        stimulus.live.fixColor = 0;
        stimulus.live.changing = task.thistrial.change;
        stimulus.live.faces = 1;
    case stimulus.seg.resp
        stimulus.live.fixColor = 1;
        stimulus.live.changing = 0;
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
upCues(task,stimulus);
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
        % down up code
% % % %         if changeTime < time / 2
% % % %             % we are still dropping
% % % %             % current frame is:
% % % %             frame = round(changeTime*2/(time/(rNum-1)))+1;
% % % %         else
% % % %             frame = (rNum*2-1)-round(changeTime*2/(time/(rNum-1)));
% % % %         end
        % down only code
        frame = round(changeTime/(time/(rNum-1)))+1;
        if frame < 1, frame = 1; end
        if frame > rNum, frame = rNum; end
        curimage = stimulus.changeTex{frame};
%         disp(sprintf('frame %i time %03.f',frame,changeTime));
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
        stimulus.staircase{task.thistrial.cues,task.thistrial.contrastList(task.thistrial.target)} = ...
            doStaircase('update',stimulus.staircase{task.thistrial.cues,task.thistrial.contrastList(task.thistrial.target)},task.thistrial.correct);
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

% Increment block (it starts at 0)
stimulus.blocks.counter = stimulus.blocks.counter + 1;

% clear screen
mglClearScreen(stimulus.colors.rmed/255);
mglTextDraw('Run starting... please wait.',[0 0]);
mglFlush

myscreen.flushMode = 1;

% let the user know
disp(sprintf('(fade) Starting block number: %i',stimulus.blocks.counter));

stimulus.blocks.curTrialPerm = stimulus.blocks.permutationMatrix(randperm(size(stimulus.blocks.permutationMatrix,1)),:);
stimulus.blocks.curContrastList = stimulus.blocks.contrastListOptions(stimulus.blocks.counter,:);
stimulus.blocks.curMaxContrast = 1;
% stimulus.blocks.curMaxContrast = max(stimulus.pedestals.contrast(stimulus.blocks.curContrastList));

% setGammaTableForMaxContrast(stimulus.blocks.curMaxContrast);

% now that we know the current contrasts values we can build the flyTex
% values for this current trial
stimulus.blocks.curGenderList = randi(2,4);
stimulus.blocks.curImageList = [-1 -1 -1 -1];
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

highCon = stimulus.pedestals.contrast(task.thistrial.contrastList(task.thistrial.changeTarget));
lowCon = highCon - task.thistrial.deltaPed;
if lowCon <= 0
    disp('(fade) Required contrast range dropped below zero, truncating');
    lowCon = .01;
end

range = fliplr(lowCon : (highCon-lowCon)/(stimulus.steps-1) : highCon);

if isfield(stimulus,'changeTex')
    for i = 1:length(stimulus.changeTex)
        mglDeleteTexture(stimulus.changeTex{i});
    end
else
    stimulus.changeTex = cell(1,length(range));
end

for i = 1:length(range)
    % get the image
    img = stimulus.images{stimulus.blocks.curGenderList(task.thistrial.changeTarget)}(:,:,stimulus.blocks.curImageList(task.thistrial.changeTarget));
    % set the contrast
    curCon = range(i);
    % scale by the ratio
    npdf = scalePdf(stimulus.basepdf,stimulus.colors.mrmin:stimulus.colors.mrmax,curCon / stimulus.blocks.curMaxContrast);
    % change the image to match the PDF
    img = (stimulus.colors.mrmax-stimulus.colors.mrmin)*histeq(img/255,npdf) + stimulus.colors.mrmin;

    stimulus.changeTex{i} = mglCreateTexture(img);
end

%% getDeltaPed

function [deltaPed, stimulus] = getDeltaPed(stimulus,tCon,cues)
[deltaPed, stimulus.staircase{cues,tCon}] = doStaircase('testValue',stimulus.staircase{cues,tCon});


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
    % set the contrast
    curCon = stimulus.pedestals.contrast(stimulus.blocks.curContrastList(i));
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

stimulus.staircase = cell(2,length(stimulus.pedestals.contrast));
stimulus.staircase{1,1} = doStaircase('init','upDown',...
        'initialThreshold',stimulus.baseThresh(1),...
        'initialStepsize',stimulus.baseThresh(1)/3,...
        'minThreshold=.001','maxThreshold=.2','stepRule','pest',...
        'nTrials=60','maxStepsize=.1','minStepsize=.001');
stimulus.staircase{2,1} = stimulus.staircase{1,1};
for i = 2:length(stimulus.pedestals.contrast)
    stimulus.staircase{1,i} = stimulus.staircase{1,1};
    stimulus.staircase{2,i} = stimulus.staircase{1,1};
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus)

try
    
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
    plot(stimulus.pedestals.contrast(1:3),plotting(1,:),'-r');
    plot(stimulus.pedestals.contrast(1:3),plotting(2,:),'--r');
    axis([stimulus.pedestals.contrast(1) stimulus.pedestals.contrast(3) 0 1]);
    hold off

catch
    disp('(fade) Figures were not generated successfully.');
end

%% checkStaircaseStop
function checkStaircaseStop(stimulus)
% Check both staircases
for cues = 1:2
    s = stimulus.staircase{cues};
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
        stimulus.staircase{cues} = doStaircase('init',s,'initialThreshold',vals{threshPos},'initialStepsize',vals{stepPos});
        STOP = 1;
    end
end

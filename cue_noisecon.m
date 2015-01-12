% cue_noisecon
%
%      usage: myscreen=cue_noisecon()
%         by: daniel birman
%       date: 11/10/14
%    purpose: contrast and noise discrimination (different blocks) cueing
%             different numbers of stimuli.
%
%        use: call cue_noisecon to initialize. The first
%             run takes significantly longer due to loading stimuli. Each
%             run is either contrast or noise. For training, tell
%             participants to perform either contrast, noise, or every
%             FIFTH run, gender categorization. Once performance stabilizes
%             (see the output from dispStaircase and dSP), tell
%             participants to perform the dual task. Make sure to set
%             dual=? during the correct runs, and let the
%             code handle the other flags.
%
%             If you want to run just the main task (e.g. in the scanner)
%             without the peripheral task, use peripheralTask=0 to disable
%             the output.
%
%      flags: peripheralTask (0/1) - Flag disables the peripheral gender
%             categorization task
%             stimFileNum (-1/#) - Load a specific stimfile from a
%             subject's folder. Defaults to the last stimfile. Warning:
%             Only the first stimfile is saved with the image file data,
%             subsequent stimfiles only contain run data.
%             dual (0/1) - Uses the dual task staircases instead of the
%             single-task staircases. Never run a participant without
%             correctly setting this flag!
%             plots (0/1) - Displays staircase plots (and estimated
%             psychophysic functions)
%             practice (def=0) - Displays on ten trials
%

function [myscreen] = cue_noisecon(varargin)

global stimulus
%% Initialize Variables

% add arguments later
widthDeg = []; heightDeg = [];
peripheralTask = [];
stimFileNum = [];
plots = [];
practice = [];
getArgs(varargin,{'widthDeg=5.5', 'heightDeg=5.5', ...
    'peripheralTask=1','stimFileNum=-1', ...
    'dual=0','plots=1','practice=0'});

stimulus.dual = dual;
stimulus.counter = 1; % This keeps track of what "run" we are on.
%% Setup Screen

screen.screenNumber = 2;
myscreen = initScreen(screen);

%% Open Old Stimfile
stimulus.initStair = 1;

if ~isempty(mglGetSID) && isdir(sprintf('~/data/cue_noisecon/%s',mglGetSID))
    % Directory exists, check for a stimefile
    files = dir(sprintf('~/data/cue_noisecon/%s/1*mat',mglGetSID));

    if length(files) >= 1
        if stimFileNum == -1
            if length(files) > 1
                warning('Multiple stimfiles found, loading last one. If you wanted a different functionality use stimFileNum=#');
            end
            fname = files(end).name;
        else
            fname = files(stimFileNum).name;
        end
        s = load(sprintf('~/data/cue_noisecon/%s/%s',mglGetSID,fname));
        stimulus.staircase = s.stimulus.staircase;
        stimulus.p.staircase = s.stimulus.p.staircase;
        stimulus.dualstaircase = s.stimulus.dualstaircase;
        stimulus.p.dualstaircase = s.stimulus.p.dualstaircase;
        stimulus.counter = s.stimulus.counter + 1;

        % load blocks too
        stimulus.blocks = s.stimulus.blocks;
        stimulus.blocks.loaded = 1;

        clear s;
        stimulus.initStair = 0;
        disp(sprintf('(noisecon) Data file: %s loaded, this is run #%i',fname,stimulus.counter));
    end
end

%% Initialize Stimulus

myscreen = initStimulus('stimulus',myscreen);

% Setup for second task
stimulus.pFlag = 0;
stimulus.pInt = 0;
stimulus.pActive = peripheralTask;

% Colors: We reserve the first few colors
stimulus.colors.nReservedPeripheral = 13;
stimulus.colors.maxPer = .6;
stimulus.colors.minPer = .4;
values = stimulus.colors.minPer:(stimulus.colors.maxPer- ...
    stimulus.colors.minPer)/(stimulus.colors.nReservedPeripheral-1):stimulus.colors.maxPer;
stimulus.colors.reservedColors = [values',values',values'];
stimulus.colors.reservedColors = [stimulus.colors.reservedColors;1 0 0;0 1 0;1 1 1;0 0 0];

stimulus.colors.nReservedColors = size(stimulus.colors.reservedColors,1);
stimulus.maxIndex = 255;
stimulus.colors.nFaceColors = stimulus.maxIndex - stimulus.colors.nReservedColors;

stimulus.colors.minFaceIndex = stimulus.maxIndex+1-stimulus.colors.nFaceColors;
stimulus.colors.maxFaceIndex = stimulus.maxIndex;

stimulus.pedestals.maxRange = (255-(stimulus.colors.nReservedColors+1))/2;


% set the reserved colors - this gives a convenient value between 0 and 1 to use the reserved colors with
for i = 1:stimulus.colors.nReservedColors
  stimulus.colors.reservedColor(i) = (i-1)/stimulus.maxIndex;
end

%% MGL Parameters

mglTextSet('Helvetica',32,stimulus.colors.reservedColor(16),...
0,0,0,0,0,0,0);
    
%% Gamma Table Initialization

% set the reserved colors
stimulus.gammaTable(1:size(stimulus.colors.reservedColors,1),1:size(stimulus.colors.reservedColors,2))=stimulus.colors.reservedColors;

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
  disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

%% Initialize Images

% set initial thresholds
stimulus.nExemplar = 5; % Number of each noise level to generate
stimulus.pedestals.contrast = [ .1 .15 .25 .4 .6 ];
baseThresh(2) = .55;
% These noise levels correspond to an SnR of 
stimulus.pedestals.noise = [ 1/(1+exp(2.5)) 1/(1+exp(1.75)) 1/(1+exp(1)) 1/(1+exp(.25)) 1/(1+exp(-.5))];
stimulus.pedestals.SnR = stimulus.pedestals.noise ./ (1-stimulus.pedestals.noise);
baseThresh(1) = .55;
stimulus.nPedestals = length(stimulus.pedestals.contrast);

% load images
stimulus.widthDeg = widthDeg;
stimulus.heightDeg = heightDeg;
stimulus.pos1 = [-3.8 +3.8 -3.8 +3.8];
stimulus.pos2 = [-3.8 -3.8 +3.8 +3.8];
categories = {'m' 'f'};
name = getenv('USER');
stimulus.p.widthPix = 162;
stimulus.p.heightPix = 193;
% pr = 450/385;
% stimulus.p.widthDeg = stimulus.widthDeg;
% stimulus.p.heightDeg = pr * stimulus.heightDeg;
stimulus.p.widthDeg = 5;
stimulus.p.heightDeg = 6;
stimulus.imageDirPer = fullfile(sprintf('/Users/%s/proj/att_awe/images/brazil_faces/',name));
stimulus.widthPix = 400;
stimulus.heightPix = 400;
stimulus.imageDirMain = fullfile(sprintf('/Users/%s/proj/att_awe/images/all_faces2/',name));
dispLoadFig = 0; keepAspectRatio = 0;

stimulus = InitStim(stimulus,categories,dispLoadFig,keepAspectRatio);

%% EYE CALIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myscreen = eyeCalibDisp(myscreen);

%% Setup Task

% This is the contrast/noise discrimination task
task{1}{1}.waitForBacktick = 1;
stimulus.seg.ITI = 1;
stimulus.seg.cue = 2;
stimulus.seg.stim1 = 3;
stimulus.seg.presp1 = 4;
stimulus.seg.stim2 = 5;
stimulus.seg.presp2 = 6;
stimulus.seg.resp = 7;
task{1}{1}.segmin = [1 1 .5 1 .5 1 1.4];
task{1}{1}.segmax = [1 1 .5 1 .5 1 1.4];
task{1}{1}.synchToVol = [0 0 0 0 0 0 0];
task{1}{1}.getResponse = [0 0 0 0 0 0 1];
task{1}{1}.randVars.calculated.blockType = nan;
% We will use the 2-4 pedestal of 1-5 options, to always have one above or
% below.
task{1}{1}.parameter.pedestal = 2:4; %% IMPORTANT %%
task{1}{1}.parameter.pedestalRandom = 2:4; %% IMPORTANT%%
task{1}{1}.parameter.dummy = 1:4; % This just makes sure the number of trials is large enough
task{1}{1}.parameter.cues = [1 4]; %% IMPORTANT %%
%%
stimulus.nPedestalOpts = length(task{1}{1}.parameter.pedestal);
task{1}{1}.random = 1;
task{1}{1}.numBlocks = 1;
if practice
    task{1}{1}.numTrials = 10;
end

task{1}{1}.randVars.calculated.interval = nan;
task{1}{1}.randVars.calculated.targetLoc = nan;
task{1}{1}.randVars.calculated.gender = nan;
task{1}{1}.randVars.calculated.genderList = nan(1,4);
task{1}{1}.randVars.calculated.noiseList = nan(4,2);
task{1}{1}.randVars.calculated.contrastList = nan(4,2);
task{1}{1}.randVars.calculated.maxContrast = nan;
task{1}{1}.randVars.calculated.deltaPed = nan;
task{1}{1}.randVars.calculated.imageNums = nan(1,4);

if isfield(stimulus,'blocks') && isfield(stimulus.blocks,'loaded')
    % We already have our blocks
    stimulus.blocks.counter = stimulus.blocks.counter + 1;
    stimulus.blocks = rmfield(stimulus.blocks,'loaded'); % remove the load field, otherwise it gets saved across runs
else
    % This is the first run, build up the blocks
    % Task switching
    types = [1 2];
    stimulus.blocks.fullBlocks = types(randperm(2));%repmat(types,1,task{1}{1}.numBlocks/2);
    stimulus.blocks.blockTypes = {'Noise','Contrast'};
    stimulus.blocks.counter = 0;
end
stimulus.blocks.curBlock = stimulus.blocks.fullBlocks(mod(stimulus.blocks.counter,2)+1);
stimulus.blocks.blockList(stimulus.blocks.counter+1) = stimulus.blocks.curBlock;

%% Full Setup
% Initialize task (note phase == 1)
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@getResponseCallback,@startTrialCallback,[],@startBlockCallback);
end
if peripheralTask
    [task{2}, myscreen] = cue_gender(myscreen);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init staircase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if stimulus.initStair
    stimulus.initThresh = zeros(2,2,3);

    for cond = 1:2
        for cues = 1:2
            stimulus.initThresh(cond,cues) = baseThresh(cond);
        end
    end
    disp(sprintf('(noisecon) Initializing staircases'));
    stimulus = initStaircase(stimulus);
else
    disp('(noisecon) Re-using staircase from previous run');

    % Reset staircase if necessary
    if stimulus.dual
        for blocks = 1:2
            for cue = 1:2
                for peds = 1:3
                    stimulus.dualstaircase{blocks,cue,peds} = checkStaircaseStop(stimulus.dualstaircase{blocks,cue,peds});
                end
            end
        end
        stimulus.p.dualstaircase{1} = checkStaircaseStop(stimulus.p.dualstaircase{1});
        stimulus.p.dualstaircase{2} = checkStaircaseStop(stimulus.p.dualstaircase{2});
    else
        for blocks = 1:2
            for cue = 1:2
                for peds = 1:3
                    stimulus.staircase{blocks,cue,peds} = checkStaircaseStop(stimulus.staircase{blocks,cue,peds});
                end
            end
        end
        stimulus.p.staircase = checkStaircaseStop(stimulus.p.staircase);
    end
end

%% Main Task Loop

phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    if peripheralTask
        [task{2}, myscreen] = updateTask(task{2},myscreen,1);
    end
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

mglClearScreen(stimulus.colors.reservedColor(7));
mglTextDraw('Run complete... please wait.',[0 0]);
mglFlush
disp('(noisecon) Run ending, displaying stimuli...');

if plots
    dispStaircase(stimulus);
    dispStaircaseP(stimulus);
end

% delete texture
if isfield(stimulus,'flyTex')
    for x = 1:size(stimulus.flyTex,1)
        for y = 1:size(stimulus.flyTex,2)
            mglDeleteTexture(stimulus.flyTex{x,y});
        end
    end
end
% delete texture
if isfield(stimulus.p,'tex')
    for x = 1:size(stimulus.p.tex,1)
        for y = 1:size(stimulus.p.tex,2)
            mglDeleteTexture(stimulus.p.tex{x,y});
        end
    end
end
    stimulus = rmfield(stimulus,'flyTex');
    stimulus.p = rmfield(stimulus.p,'tex');

if stimulus.counter > 1
    % Temporarily remove raw files to limit stimfile size, then replace these
    % so they can be used on the next run.
    stimbackup = stimulus;
    stimulus = rmfield(stimulus,'raw');
    stimulus = rmfield(stimulus,'images');
    stimulus.p = rmfield(stimulus.p,'raw');
    stimulus.p = rmfield(stimulus.p,'images');
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

global stimulus

stimulus.curTrial = task.trialnum;
stimulus.values = zeros(4,2);

disp(sprintf('(cue_noisecon) Starting trial %i',task.blockTrialnum));

% Set the basic missing parameters
pos_int = 1:2;
task.thistrial.interval = pos_int(randi(2));
pos_loc = 1:4;
task.thistrial.targetLoc = pos_loc(randi(4));
pos_gen = 1:2;
task.thistrial.gender = pos_gen(randi(2));

% Let's get four random images, with no repetitions
imgs = zeros(1,4);
while any(imgs==0) || length(unique(imgs))<length(imgs)
    imgs = round(rand(1,4)*stimulus.raw{1}.n);
end
task.thistrial.imageNums = imgs;

pos = [1 2 3 4];
% Let's set the distractor image pedestals (the other three heights)  
otherPedestals = [randi(stimulus.nPedestals) randi(stimulus.nPedestals) randi(stimulus.nPedestals)];
otherPedestals = otherPedestals(randperm(length(otherPedestals)));

pedestals(pos==task.thistrial.targetLoc) = task.thistrial.pedestal;
pedestals(pos~=task.thistrial.targetLoc) = otherPedestals;

% We also need the randoms, if we are in a noise block then noise will get
% the pedestals and contrast will get the randoms.
rands = [randi(stimulus.nPedestals) randi(stimulus.nPedestals) randi(stimulus.nPedestals)];
rands = rands(randperm(length(rands)));

randoms(pos==task.thistrial.targetLoc) = task.thistrial.pedestalRandom;
randoms(pos~=task.thistrial.targetLoc) = rands;

% Let's set the image genders

task.thistrial.genderList = [0 0 0 0];
task.thistrial.genderList(task.thistrial.targetLoc) = task.thistrial.gender;
if task.thistrial.gender == 1
    gens = [1 2 2];
else
    gens = [1 1 2];
end
task.thistrial.genderList(task.thistrial.genderList==0) = gens(randperm(3));

% Get details
blocks = stimulus.blocks.curBlock;
cue = find(task.thistrial.cues==[1 4]);
peds = task.thistrial.pedestal-1;

% Get Delta
[task.thistrial.deltaPed, stimulus] = getDeltaPed(stimulus,blocks,cue,peds);

    
% Get maxContrast
task.thistrial.maxContrast = getMaxContrast(task,stimulus,task.thistrial.pedestal,pedestals,randoms);
if task.thistrial.maxContrast > 1, task.thistrial.maxContrast = 1; end % Make sure we don't exceed 1
setGammaTableForMaxContrast(task.thistrial.maxContrast);

% Set up each of the image textures into stimulus.flyTex
for imagePos = 1:4
    p_mask = rand(1,length(stimulus.raw{1}.halfFourier{1}.phase))*2*pi;
    for int = 1:2
        if imagePos==task.thistrial.targetLoc || int==1
            [task, stimulus.flyTex{imagePos,int}] = convertToTex(imagePos,task,imagePos==task.thistrial.targetLoc,int,pedestals(imagePos),randoms(imagePos),p_mask);
        else
            % We can save a tiny bit of time by copying the image from the
            % other interval if this isn't the target image.
            stimulus.flyTex{imagePos,int} = stimulus.flyTex{imagePos,1};
        end
    end
end

%% CONVERT TO TEX
function [task, tex] = convertToTex(imgNum,task,isTarget,int,ped,pedR,p_mask)
global stimulus

img = stimulus.images{task.thistrial.genderList(imgNum)}(:,:,task.thistrial.imageNums(imgNum));

add = 0;
if isTarget && int==task.thistrial.interval
    % Get the pedestal delta note, this value depends on blocktype and is
    % specific to contrast or noise.
    add = task.thistrial.deltaPed;
end
% First we need to know whether we are adjusting noise or contrast
% randomly.
if stimulus.blocks.curBlock == 1
    %%%% PEDESTAL = NOISE %%%%
    curNoi = stimulus.pedestals.noise(ped) + add;
    curCon = stimulus.pedestals.contrast(pedR);
else
    %%%% PEDESTAL = CONTRAST %%%%
    curNoi = stimulus.pedestals.noise(pedR);
    curCon = stimulus.pedestals.contrast(ped) + add;
end
[curNoi, curCon] = checkValues(curNoi,curCon);

task.thistrial.noiseList(imgNum,int) = curNoi;
task.thistrial.contrastList(imgNum,int) = curCon;

img_f = getHalfFourier(img);

% Add NOISE
img = pinkNoise(img,img_f,curNoi,p_mask);

% Add CONTRAST
% we know maxContrast, get the range values
rmed = stimulus.colors.nReservedColors + 1 + stimulus.pedestals.maxRange;
mrmax = rmed + stimulus.pedestals.maxRange;
mrmin = rmed - stimulus.pedestals.maxRange;
% build the normalized PDF
npdf = normpdf(mrmin:mrmax,rmed,50);
% scale by the ratio
npdf = scalePdf(npdf,mrmin:mrmax,curCon / task.thistrial.maxContrast);
% change the image to match the PDF
img = (mrmax-mrmin)*histeq(img/255,npdf) + mrmin;

% imageSaver(img,imgNum,curCon,curNoi,task.thistrial.maxContrast);

tex = mglCreateTexture(img);

%% IMG SAVE

function imageSaver(img,num,c,n,m)
warning('IMAGE SAVER ENABLED!!!');
if m==.8
    if ~isdir(sprintf('~/data/imageExamples/%s',num2str(num)))
        mkdir(sprintf('~/data/imageExamples/%s',num2str(num)));
    end
    imwrite(flipud(img/255),sprintf('~/data/imageExamples/%s/image%i_c%0.2f_n%0.2f_m%0.2f.tif',num2str(num),c,n,m),'tif');
end

%% getMaxContrast

function maxC = getMaxContrast(task,stimulus,cPed,peds,rands)

if stimulus.blocks.curBlock == 1
    % Pedestal = noise, so the max contrast will be the maximum value of
    % the randoms
    maxC = stimulus.pedestals.contrast(max(rands));
else
    % Pedestal = contrast, so the max contrast will be the maximum value of
    % either the top contrast
    maxC = max(stimulus.pedestals.contrast(cPed)+task.thistrial.deltaPed,max(peds));
end

function [cN, cC] = checkValues(cN,cC)
if cC > 1
    warning('Max contrast exceeded 1. Thresholding');
    cC = 1;
end
if cN > 1
    warning('Max noise exceeded 1. Thresholding');
    cN = 1;
end
%% getDeltaPed

function [deltaPed, stimulus] = getDeltaPed(stimulus,condition,cue,p)
if stimulus.dual
    [deltaPed, stimulus.dualstaircase{condition,cue,p}] = doStaircase('testValue',stimulus.dualstaircase{condition,cue,p});
else
    [deltaPed, stimulus.staircase{condition,cue,p}] = doStaircase('testValue',stimulus.staircase{condition,cue,p});
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Runs at the start of each Segment %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = startSegmentCallback(task, myscreen)

global stimulus
myscreen.flushMode = 0;
if any(task.thistrial.thisseg == [stimulus.seg.stim1 stimulus.seg.stim2])
    stimulus.pFlag = 1;
    stimulus.pInt = (task.thistrial.thisseg - 1) / 2;
elseif any(task.thistrial.thisseg == [stimulus.seg.resp stimulus.seg.ITI stimulus.seg.cue stimulus.seg.presp1 stimulus.seg.presp2])
    stimulus.pFlag = 2;
else
    stimulus.pFlag = 0;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Refreshes the Screen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task, myscreen] = screenUpdateCallback(task, myscreen)
global stimulus

if myscreen.flushMode == 0
    mglClearScreen(stimulus.colors.reservedColor(7));
    stimulus.fixColor = stimulus.colors.reservedColor(17);
    % set the fixation color
    if any(task.thistrial.thisseg == [stimulus.seg.cue stimulus.seg.presp1 stimulus.seg.presp2])
        % Cue segment or the peripheral response segments
        upCues(task,stimulus);
        upFix(stimulus);
    elseif task.thistrial.thisseg == stimulus.seg.resp
        % Resp -- uses bright fixation cross
        upCues(task,stimulus);
        stimulus.fixColor = stimulus.colors.reservedColor(16);
        upFix(stimulus);
    elseif any(task.thistrial.thisseg == [stimulus.seg.stim1 stimulus.seg.stim2])
        % Either of the stimulus segments
        upCues(task,stimulus);
        upFaces(stimulus,task);
        % This turned the fixation color white, but it seems unnecessarily
        % distracting.
        upFix(stimulus);
    else
        upFix(stimulus);
    end
end

%%
function upFix(stimulus)

mglFixationCross(1,1,stimulus.fixColor);

%%
function upCues(task,stimulus)

ang = [180 0 180 0];
if task.thistrial.cues == 1 || task.thistrial.thisseg == stimulus.seg.resp
    % We can't just display all of the lines, we just want one line
    usePos1 = stimulus.pos1(task.thistrial.targetLoc);
    usePos2 = stimulus.pos2(task.thistrial.targetLoc);
    ang = ang(task.thistrial.targetLoc);
else
    usePos1 = stimulus.pos1;
    usePos2 = stimulus.pos2;
end
ang = d2r(ang);
mglLines2(cos(ang + atan(usePos1./usePos2))*.1,sin(ang + atan(usePos1./usePos2))*.1, ...
    cos(ang + atan(usePos1./usePos2))*.75,sin(ang + atan(usePos1./usePos2))*.75,1,stimulus.fixColor);

function upFaces(stimulus,task)

for imagePos = 1:4
    % Get the image from the textures built on the fly
    curimage = stimulus.flyTex{imagePos,(task.thistrial.thisseg-1)/2};
    % Move to the screen
    mglBltTexture(curimage,[stimulus.pos1(imagePos) stimulus.pos2(imagePos) stimulus.widthDeg stimulus.heightDeg]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Called When a Response Occurs %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = getResponseCallback(task, myscreen)

global stimulus
mglClearScreen(stimulus.colors.reservedColor(7));
% stimulus.imagest = num2str(task.thistrial.thisseg);
% upText(stimulus);

%%%% NOTE: PARTICIPANTS ALWAYS LOOK FOR THE MORE VISIBLE ITEM, i.e. HIGHER NOISE OR CONTRAST %%%%

block = stimulus.blocks.curBlock;
cue = find(task.thistrial.cues==[1 4]);
pedPos = task.thistrial.pedestal-1;    

if any(task.thistrial.whichButton == [1 2])
    if task.thistrial.gotResponse == 0
        whichInterval = find(task.thistrial.interval == [1 2]);
        if (task.thistrial.whichButton == whichInterval)
            correctIncorrect = 'correct';
            stimulus.fixColor = stimulus.colors.reservedColor(15);
            if stimulus.dual
                stimulus.dualstaircase{block,cue,pedPos} = ...
                    doStaircase('update',stimulus.dualstaircase{block,cue,pedPos},1);
            else
                stimulus.staircase{block,cue,pedPos} = ...
                    doStaircase('update',stimulus.staircase{block,cue,pedPos},1);
            end
        else
            correctIncorrect = 'incorrect';
            stimulus.fixColor = stimulus.colors.reservedColor(14);  
            if stimulus.dual
                stimulus.dualstaircase{block,cue,pedPos} = ...
                    doStaircase('update',stimulus.dualstaircase{block,cue,pedPos},0);
            else
                stimulus.staircase{block,cue,pedPos} = ...
                    doStaircase('update',stimulus.staircase{block,cue,pedPos},0);
            end
        end
        disp(sprintf('(noisecon) Response %s',correctIncorrect));
    else
        disp(sprintf('Subject responded multiple times: %i',task.thistrial.gotResponse+1));
    end
end
myscreen.flushMode = 1;
upFix(stimulus);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Block Call Back %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [task, myscreen] = startBlockCallback(task, myscreen)

global stimulus

mglClearScreen(stimulus.colors.reservedColor(7));

stimulus.blockList(stimulus.counter) = stimulus.blocks.curBlock;
stimulus.dualList(stimulus.counter) = stimulus.dual;

myscreen.flushMode = 1;
mglTextDraw(stimulus.blocks.blockTypes{stimulus.blocks.curBlock},[0,0]);
disp(sprintf('(noisecon) Block type: %s',stimulus.blocks.blockTypes{stimulus.blocks.curBlock}));

%% scalePDF

function [scaled, X2] = scalePdf(unscaled,X,factor)
if factor > 1
    error('Scale not designed for scale > 1');
end

X2 = min(X):factor*(X(2)-X(1)):max(X);
pad = floor((length(X2)-length(X))/2);
scaled = [zeros(1,pad) unscaled zeros(1,pad)];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pinkNoise %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function image = pinkNoise(img,img_f,K,p_mask)

% Get the image range
ma = max(img(:));
mi = min(img(:));

if K == 1
    image = img;
    return
end

% Generate white noise
nimg_f = img_f;
nimg_f.phase = p_mask;
noise = reconstructFromHalfFourier(nimg_f);

L0 = mean2(img);
image = L0 + sqrt(K) * (img - L0) + sqrt(1-K) * (noise - L0);

% Make sure we are inside range
image(image>ma) = ma;
image(image<mi) = mi;


%%%%%%%%%%%%%%%%%%%%%%%%
%    initStaircase     %
%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = initStaircase(stimulus)

stimulus.staircase = cell(2,2,stimulus.nPedestalOpts);
for condition = 1:2 % noise / contrast
    for cues = 1:2
        stimulus.staircase{condition,cues,1} = doStaircase('init','upDown', ...
            'initialThreshold',stimulus.initThresh(condition,cues), ...
            'initialStepsize',stimulus.initThresh(condition,cues)/3, ...
            'minThreshold=0.001','maxThreshold=.55','stepRule','levitt', ...
            'nTrials=50');
        stimulus.staircase{condition,cues,2} = stimulus.staircase{condition,cues,1};
        stimulus.staircase{condition,cues,3} = stimulus.staircase{condition,cues,1};
        stimulus.dualstaircase{condition,cues,1} = stimulus.staircase{condition,cues,1};
        stimulus.dualstaircase{condition,cues,2} = stimulus.staircase{condition,cues,1};
        stimulus.dualstaircase{condition,cues,3} = stimulus.staircase{condition,cues,1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircase    %
%%%%%%%%%%%%%%%%%%%%%%%
function dispStaircase(stimulus)

try
    if stimulus.dual
        stairtype = 'dualstaircase';
    else
        stairtype = 'staircase';
    end

    plotting = zeros(2,3);
    if stimulus.blocks.curBlock == 1
        % noise
        if isfield(stimulus.pedestals,'SnR')
            typeP = 'SnR';
        else
            typeP = 'noise';
        end
        num = 1;
    else
        typeP = 'contrast';
        num = 2;
    end
    
    figure % this is the 'staircase' figure
    title(sprintf('%s, Staircase plot (R->G->B high)',stimulus.blocks.blockTypes{stimulus.blocks.curBlock}));
    hold on
    drawing = {'-r' '-g' '-b'
                '--r' '--g' '--b'};
    for cues = 1:2
        for ped = 1:3
            try
                testV = [];
                for i = 1:length(stimulus.(stairtype){num,cues,ped})
                    testV = [testV stimulus.(stairtype){num,cues,ped}(i).testValues];
                end
                plot(testV,drawing{cues,ped});
            catch
            end
            try
                out = doStaircase('threshold',stimulus.(stairtype){num,cues,ped},'type','weibull'); % noise, 1 cue, lowest
                plotting(cues,ped) = out.threshold;
            catch
                plotting(cues,ped) = -1;
            end
        end
    end
    hold off
    figure
    hold on
    title(sprintf('%s, R->G->B High',stimulus.blocks.blockTypes{stimulus.blocks.curBlock}));
    plot(stimulus.pedestals.(typeP)(2:4),plotting(1,:),'-r');
    plot(stimulus.pedestals.(typeP)(2:4),plotting(2,:),'--r');
    axis([stimulus.pedestals.(typeP)(1) stimulus.pedestals.(typeP)(5) 0 1]);
    hold off

catch
    disp('(noisecon) Figures were not generated successfully.');
end

%%%%%%%%%%%%%%%%%%%%%%%
%    dispStaircaseP    %
%%%%%%%%%%%%%%%%%%%%%%
function dispStaircaseP(stimulus)

try
    if stimulus.dual
        doStaircase('threshold',stimulus.p.dualstaircase{stimulus.blocks.curBlock},'dispFig',1,'type','weibull');
    %     title('Gender Task -- estimated Threshold');
    else
        doStaircase('threshold',stimulus.p.staircase,'dispFig',1,'type','weibull');
    %     title('Gender Task (DUAL) -- estimated Threshold');
    end
catch
    disp('(gender) Figures were not generated successfully.');
end

%% checkStaircaseStop
function [s] = checkStaircaseStop(s)
if doStaircase('stop',s)
    % this is a bit of a pain... you can't pass an initialThreshold
    % argument do doStaircase('init',s, ...), it ignores everything and
    % resets using the calculated threshold. Because you can't override it
    [args, vals, ~] = getArgs(s(1).initArgs);
    if length(args) > 6
        error('Not designed for initArgs > 6, check your code!');
    end
    s(end+1) = doStaircase('init','updown',args{1},vals{1},args{2},vals{2},args{3},vals{3},args{4},vals{4},args{5},vals{5},args{6},vals{6});
end
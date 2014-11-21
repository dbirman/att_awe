function myscreen = stimulus_tests()


screen.screenNumber = 2;
myscreen = initScreen(screen);

global stimulus
myscreen = initStimulus('stimulus',myscreen);
widthPix=250;
heightPix=250;
widthDeg=12;heightDeg=12;

stimulus.colors.nReservedPeripheral = 40;
values = 0:1/stimulus.colors.nReservedPeripheral:1;
stimulus.colors.reservedColors = [values',values',values'];
stimulus.colors.reservedColors = [stimulus.colors.reservedColors;1 1 1;1 0 0;0 1 0];

stimulus.colors.nReservedColors = size(stimulus.colors.reservedColors,1);
stimulus.maxIndex = 255;
stimulus.colors.nFaceColors = stimulus.maxIndex - stimulus.colors.nReservedColors;

stimulus.colors.minFaceIndex = stimulus.maxIndex-stimulus.colors.nFaceColors+1;
stimulus.colors.maxFaceIndex = stimulus.maxIndex;

%% Initialize Images

texp = 15;

% set initial thresholds
% stimulus.pedestals.contrast = .1:.2:.9;
stimulus.pedestals.noise = fliplr(texp.^(0:.2:1)/texp-1/texp);
stimulus.pedestals.contrast = .2:.2:1;

% stimulus range
stimulus.pedestals.maxRange = (255-(stimulus.colors.nReservedColors+1))/2;

stimulus.gT = 0;

% load images
stimulus.widthPix = widthPix;
stimulus.heightPix = heightPix;
stimulus.widthDeg = widthDeg;
stimulus.heightDeg = heightDeg;
categories = {'m','f'};
% objlocDir = which('objloc');objlocDir = fileparts(objlocDir);
% imageDir = fullfile(objlocDir,'images/ObjLocImages');
imageDir = fullfile('/Users/dan/proj/att_awe/images/af_small/');
dispLoadFig = 0; keepAspectRatio = 0;

% set the reserved colors
stimulus.gammaTable(1:size(stimulus.colors.reservedColors,1),1:size(stimulus.colors.reservedColors,2))=stimulus.colors.reservedColors;

% set the reserved colors - this gives a convenient value between 0 and 1 to use the reserved colors with
for i = 1:stimulus.colors.nReservedColors
  stimulus.colors.reservedColor(i) = (i-1)/stimulus.maxIndex;
end

saveFile = fullfile('/Users/dan/proj/att_awe/s_t_sf.mat');
stimulus = InitStim(stimulus,myscreen,categories,imageDir,dispLoadFig,keepAspectRatio,saveFile);

%% Gamma Table Initialization

% get gamma table
if ~isfield(myscreen,'gammaTable')
  stimulus.linearizedGammaTable = mglGetGammaTable;
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
  disp(sprintf('(cuecon:initGratings) No gamma table found in myscreen. Contrast displays like this'));
  disp(sprintf('         should be run with a valid calibration made by moncalib for this monitor.'));
  disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
end
stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;

%% Testing

% This is the contrast/noise discrimination task
task{1}{1}.waitForBacktick = 1;
task{1}{1}.parameter.noise = 1:length(stimulus.pedestals.noise);
task{1}{1}.parameter.contrast = 1:length(stimulus.pedestals.contrast);
task{1}{1}.parameter.exemplar = 1:10;
task{1}{1}.segmin = [.05 1];
task{1}{1}.segmax = [.05 1];
% task{1}{1}.parameter.exemplar = 1:10;
task{1}{1}.parameter.gender = [1 2];
task{1}{1}.parameter.number = 1:length(stimulus.tex{1,1,1});
task{1}{1}.random = 0;


for phaseNum = 1:length(task{1})
    [task{1}{phaseNum}, myscreen] = initTask(task{1}{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,[],@startTrialCallback,[],@startBlock);
end


phaseNum = 1;
% Again, only one phase.
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
    % update the task
    [task{1}, myscreen, phaseNum] = updateTask(task{1},myscreen,phaseNum);
    % flip screen
    myscreen = tickScreen(myscreen,task);
end

clear stimulus

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%

function [task myscreen] = startTrialCallback(task, myscreen)
global stimulus
myscreen.flushMode = 0;
% Now setup the contrast profiles
stimulus.gT = rand<.5; % are we using the gamma table?

cC = stimulus.pedestals.contrast(task.thistrial.contrast);
if stimulus.gT
    setGammaTableForMaxContrast(cC);
else
    setGammaTableForMaxContrast(1);
end
%%

function [task myscreen] = startBlock(task, myscreen)

global stimulus

stimulus.p_mask = rand(1,31000)*2*pi;


%%

function [task myscreen] = startSegmentCallback(task, myscreen)
global stimulus

if task.thistrial.thisseg==1
    mglClearScreen(stimulus.colors.reservedColor(15));
    return
end

mglClearScreen(stimulus.colors.reservedColor(15));

try
    mglDeleteTexture(stimulus.curTex);
end

% Get the current image number
stimulus.imageN = task.thistrial.number;
% Get the logical array for noise position
cN = stimulus.pedestals.noise(task.thistrial.noise);
% Set the gamma table
cC = stimulus.pedestals.contrast(task.thistrial.contrast);

% note: stimulus.tex{GENDER,RANGE}(:,:,#)
curImage = stimulus.tex{task.thistrial.gender}(:,:,stimulus.imageN);

% Before converting to a texture, do the noise manipulation
curImage_f = getHalfFourier(curImage);
p_mask = stimulus.p_mask(:);
p_mask = p_mask(1:length(curImage_f.phase))';
% Add noise
curImage = pinkNoise(curImage,curImage_f,cN,p_mask);

rmed = stimulus.colors.nReservedColors + 1 + stimulus.pedestals.maxRange;
mrmax = rmed + stimulus.pedestals.maxRange;
mrmin = rmed - stimulus.pedestals.maxRange;

if stimulus.gT
    % Give the image the full possible range
    npdf = normpdf(mrmin:mrmax,rmed,50);
else
    % Mess with the image
    npdf = normpdf(mrmin:mrmax,rmed,50);
    % This time scale the PDF
    npdf = scalePdf(npdf,mrmin:mrmax,cC);
end
% Build the final image
curImage2 = (mrmax-mrmin)*histeq(curImage/255,npdf) + mrmin;

% Send the texture to be BLT
stimulus.curTex = mglCreateTexture(curImage2);

% Add this face to its position
mglBltTexture(stimulus.curTex,[0 0 stimulus.widthDeg stimulus.heightDeg]);
% Mask this, based on the mask defined in initStim
% % % % % % mglBltTexture(stimulus.mask,[0 0 stimulus.widthDeg stimulus.heightDeg]);
myscreen.flushMode = 1;
mglTextDraw(sprintf('Noise: %0.2f, Contrast: %0.2f',cN,cC),[0,-7]);
mglTextDraw(sprintf('Adjusting Gamma: %i',stimulus.gT),[0,-8]);

function [task myscreen] = screenUpdateCallback(task, myscreen)

% do nothing

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
function stimulus = NC_InitStimulus(stimulus,myscreen,categories,imageDir,dispFig,keepAspectRatio)

if ~isfield(stimulus,'cuenoise'),stimulus.imagesLoaded = 0;end
stimulus.cuenoise = 1;
stimulus.thisTex = [];

if ~isfield(stimulus,'nReservedColors')
    stimulus.nReservedColors = size(stimulus.colors.reservedColors,1);
end

%% Image Loading

% make sure widht and height are odd
if iseven(stimulus.widthPix), stimulus.widthPix = stimulus.widthPix-1;end
if iseven(stimulus.heightPix), stimulus.heightPix = stimulus.heightPix-1;end

% check whether images are loaded
averageN = 0;
if ~isfield(stimulus,'imagesLoaded') || (~stimulus.imagesLoaded) || ~isequal(stimulus.categories,categories) || ~isequal(stimulus.imageDir,imageDir) || dispFig
  stimulus.nCategories = length(categories);

  % keep the averageMag and averageDc so that we can normalize images
  stimulus.averageMag = 0;
  stimulus.averageDC = 0;

  for i = 1:stimulus.nCategories
    if ~any(strcmp(categories{i},{'scramble','blank','gray'}))
      % load images
      stimulus.raw{i} = loadNormalizedImages(fullfile(imageDir,categories{i}),'width',stimulus.widthPix,'height',stimulus.heightPix,'dispFig',dispFig,'keepAspectRatio',keepAspectRatio);

      % make sure we opened ok
      if isempty(stimulus.raw{i})
	disp(sprintf('(objloc) Could not load images; %s',categories{i}));
	keyboard
      end
      % get average mag
      stimulus.averageMag = stimulus.averageMag + stimulus.raw{i}.averageMag;
      stimulus.averageDC = stimulus.averageDC +stimulus.raw{i}.averageDC;
      averageN = averageN + 1;
    else
      stimulus.raw{i}.n = 0;
    end
  end
  % compute average
  stimulus.averageMag = stimulus.averageMag/averageN;
  stimulus.averageDC = stimulus.averageDC/averageN;
  % and set that we have loaded
  stimulus.imageDir = imageDir;
  stimulus.imagesLoaded = 1;

else
  disp(sprintf('(objloc) Stimulus already initialized'));
end

% count how many images
nImages = [];
for i = 1:stimulus.nCategories
  if stimulus.raw{i}.n > 0
    % keep number of images so that we can make scrambles with equal number of images
    nImages(end+1) = stimulus.raw{i}.n;
  end
end

%% Contrast Prep
% set maximum color index (for 24 bit color we have 8 bits per channel, so 255)
% maxIndex = 255;
% 
% % get gamma table
% if ~isfield(myscreen,'gammaTable')
%   stimulus.linearizedGammaTable = mglGetGammaTable;
%   disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')); %#ok<*DSPS>
%   disp(sprintf('(taskTemplateContrast10bit:initGratings) No gamma table found in myscreen. Contrast'));
%   disp(sprintf('         displays like this should be run with a valid calibration made by moncalib'));
%   disp(sprintf('         for this monitor.'));
%   disp(sprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'));
% end
% stimulus.linearizedGammaTable = myscreen.initScreenGammaTable;
% 
% % number of colors possible for gratings, make sure that we 
% % have an odd number
% stimulus.colors.nGratingColors = maxIndex+1-stimulus.colors.nReservedColors;
% if iseven(stimulus.colors.nGratingColors)
%   stimulus.colors.nGratingColors = stimulus.colors.nGratingColors-1;
% end


%% Build Texture array

% note: stimulus.tex{PHASE,CONTRAST}(#)
stimulus.LUT = [];

for cat = 1:length(stimulus.raw)
    for N = 1:length(stimulus.pedestals.noise)
        for C = 1:length(stimulus.pedestals.contrast)
            for imgN = 1:stimulus.raw{cat}.n
                
                cN = stimulus.pedestals.noise(N);
                cC = stimulus.pedestals.contrast(C);
                
                thisImage = stimulus.raw{cat}.halfFourier{imgN};
                thisImage.dc = stimulus.averageDC;
                thisImage.mag = stimulus.averageMag;
                
                %% MASK
                img_mask = stimulus.raw{cat}.im(:,:,imgN)==0;
                
                %% NOISE
%                 thisImage = phaseScramble(thisImage,img_mask,cN);
                thisImage = pinkNoise(thisImage,img_mask,cN);
                
                %% CONTRAST
%                 thisImage.mag = thisImage.mag * cC;
                
                %% Add everything to stimulus.tex{Noise,Contrast}(img)
%                 thisImage = reconstructFromHalfFourier(thisImage);
                
                stimulus.tex{cat,N,C}(imgN) = mglCreateTexture(contrastNormalize(flipud(thisImage)));
                disppercent(calcPercentDone(cat,length(stimulus.raw),N,length(stimulus.pedestals.noise),C,length(stimulus.pedestals.contrast),imgN,stimulus.raw{cat}.n));
                % Also build the lookup table...
                stimulus.LUT(N,C,:) = [cN cC];
            end
        end
    end
end

disppercent(inf);

stimulus.trialStart = mglGetSecs;
stimulus.categories = categories;

function image = pinkNoise(img_f,mask,K)

img = reconstructFromHalfFourier(img_f);

% Generate white noise
noise = 255*spatialPattern(img_f.originalDims,-1);

L0 = mean2(img + noise) / 2; %% TODO: This should really be avoided, since all the images will have different luminance based on the noise.

image = zeros(img_f.originalDims);
output = L0 + sqrt(K) * (img - L0) + sqrt(1-K) * (noise - L0);
image(mask==0) = output(mask==0);
image(mask==1) = 76.5;

function image = phaseScramble(img_f,mask2,cN)

image = reconstructFromHalfFourier(img_f);
image = im_wmpscramble(image,cN);
image(mask2==1) = 76.5;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %    contrastNormalize    %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
function im = contrastNormalize(im)
% 
% % new version (July 2011) to keep luminance and rms contrast
% % constant across images. Scaling factor 0.7 chosen for a specific set
% % of stimuli (the set of images that was in grustim/images/ObjLocImages/
% % in June 2011)
% 
% %  im = im*0.7;
% %  im(im>255) = 255;
% %  im(im<0) = 0;

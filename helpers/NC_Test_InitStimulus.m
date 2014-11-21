function stimulus = NC_Test_InitStimulus(stimulus,myscreen,categories,imageDir,dispFig,keepAspectRatio)

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
      stimulus.raw{i} = loadNormalizedEqImages(fullfile(imageDir,categories{i}),'width',stimulus.widthPix,'height',stimulus.heightPix,'dispFig',dispFig,'keepAspectRatio',keepAspectRatio);

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

%% Build Texture array

% note: stimulus.tex{PHASE,CONTRAST}(#)
stimulus.LUT = [];

for cat = 1:length(stimulus.raw)
    for N = 1:length(stimulus.pedestals.noise)
        for C = 1:length(stimulus.pedestals.contrast)
            for imgN = 1:stimulus.raw{cat}.n
                
                cN = stimulus.pedestals.noise(N);
                cC = stimulus.pedestals.contrast(C);
                
                %% MASK
                img_mask = stimulus.raw{cat}.im(:,:,imgN)==0;
                
                %% Get IMAGE
                
                thisImage = stimulus.raw{cat}.halfFourier{imgN};
                image2 = reconstructFromHalfFourier(thisImage);
                thisImage.dc = stimulus.averageDC;
                
                %% Smooth the average magnitude CDF
                image1 = reconstructFromHalfFourier(thisImage);
                thisImage.mag = stimulus.averageMag;
                
                image = reconstructFromHalfFourier(thisImage);
%                 [avg_CDF x] = hist(sqrt(stimulus.averageMag),64);
%                 avg_CDF = cumsum(avg_CDF);
%                 p = sigm_fit(x,avg_CDF,[0, 80000, NaN, NaN]);
%                  
                
                %% NOISE
                for exemplar = 1:5
%                     exemImg = pinkNoise(thisImage,img_mask,cN);
                    exemImg = image;
                    stimulus.tex{cat,N,C,exemplar}(imgN) = mglCreateTexture(contrastNormalize(flipud(exemImg)));
                end
                
                %% CONTRAST
%                 thisImage.mag = thisImage.mag * cC;
                
                %% Add everything to stimulus.tex{Noise,Contrast}(img)
                
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
noise = 255*spatialPattern(img_f.originalDims,-.5);

im_masked = img(mask==0); % actual image;
imCDF = cumsum(histc(im_masked,0:5:255));
no_masked = noise(mask==0);
noCDF = cumsum(histc(no_masked,0:5:255));

LU = interp1(imCDF,0:255,noCDF);

[ref_counts, ~] = hist(im_masked,64);

no_masked_edit = histeq(no_masked,ref_counts);


L0 = mean2(img + noise) / 2;

image = zeros(img_f.originalDims);
output = L0 + sqrt(K) * (img - L0) + sqrt(1-K) * (noise - L0);
image(mask==0) = output(mask==0);
image(mask==1) = 0;

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
 im = im*0.7;
 im(im>255) = 255;
 im(im<0) = 0;

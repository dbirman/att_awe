function stimulus = InitStimFade(stimulus,categories,dispFig,keepAspectRatio)

%% Image Loading

% make sure widht and height are odd
if iseven(stimulus.widthPix), stimulus.widthPix = stimulus.widthPix-1;end
if iseven(stimulus.heightPix), stimulus.heightPix = stimulus.heightPix-1;end

% Normalize the 10-level image by the equalized histogram
rmed = .5*255;
mrmax = 255;
mrmin = 0;
% build the normalized PDF
npdf = normpdf(mrmin:mrmax,rmed,80);
npdf = npdf / sum(npdf);

% check whether images are loaded
averageN = 0;
if ~isfield(stimulus,'imagesLoaded') || (~stimulus.imagesLoaded) || ~isequal(stimulus.categories,categories) || dispFig
  stimulus.nCategories = length(categories);

  % keep the averageMag and averageDc so that we can normalize images
  stimulus.averageMag = 0;
  stimulus.averageDC = 0;

  for i = 1:stimulus.nCategories
    if ~any(strcmp(categories{i},{'scramble','blank','gray'}))
      % load images
      stimulus.raw{i} = loadNormEqImages2(fullfile(stimulus.imageDirMain,categories{i}),'width',stimulus.widthPix,'height',stimulus.heightPix,'dispFig',dispFig,'keepAspectRatio',keepAspectRatio);

      % make sure we opened ok
      if isempty(stimulus.raw{i})
        disp(sprintf('(initstim) Could not load images; %s',categories{i}));
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
  stimulus.imageDir = stimulus.imageDirMain;
  
  disp(sprintf('(initstim) Building main images... ___'));
  %% MAIN IMAGES
  for cat = 1:length(stimulus.raw)
      for imgN = 1:stimulus.raw{cat}.n
          thisImage = stimulus.raw{cat}.halfFourier{imgN};
          thisImage.dc = stimulus.averageDC;
          thisImage.mag = stimulus.averageMag;
          
          image = reconstructFromHalfFourier(thisImage);
          
          image = fixBoundaries(image,stimulus.pedestals.maxRange);
          stimulus.images{cat}(:,:,imgN) = flipud(image);
          
          % Display percent done
          disppercent(calcPercentDone(cat,length(stimulus.raw),imgN,stimulus.raw{cat}.n));
      end
  end
  
  stimulus.imagesLoaded = 1;
else
  disp(sprintf('(initstim) Stimulus already initialized'));
end

disppercent(inf);

stimulus.trialStart = mglGetSecs;
stimulus.categories = categories;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        pinkNoise       %
%%%%%%%%%%%%%%%%%%%%%%%%%%

function image = pinkNoise(img,img_f,K)

% Get the image range
ma = max(img(:));
mi = min(img(:));

if K == 1
    image = img;
    return
end

% Generate white noise
img_f.phase = d2r(rand(size(img_f.phase))*360-180);
noise = reconstructFromHalfFourier(img_f);

L0 = mean2(img);
image = L0 + sqrt(K) * (img - L0) + sqrt(1-K) * (noise - L0);

% Make sure we are inside range
image(image>ma) = ma;
image(image<mi) = mi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    replaceColors        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function image2 = replaceColors(image,N)
N = ceil(N);
qts = 1/N:1/N:1;
cutoffs = fliplr(quantile(image(:),qts));
image2 = image;
for i = 1:N
    image2(image<=cutoffs(i)) = N-i;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fixBoundaries        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rescale an image with any range to be within 0->255 + offset from the
% number of reserved colors.

function image = fixBoundaries(image,range)
global stimulus
if range > stimulus.pedestals.maxRange
    error('Requested range is larger than maximum range!');
end

% Get original image boundaries
ma = max(image(:));
mi = min(image(:));

% Get the new boundaries
rmed = stimulus.colors.nReservedColors + 1 + stimulus.pedestals.maxRange;
rmax = rmed + range;
rmin = rmed - range;

if rmin < stimulus.colors.nReservedColors
    error('Boundary fixing failed, image contrast range is too large');
end

% Scale image to full range
image = image * (rmax - rmin) / (ma - mi);
% Move to the image median
image = image - min(image(:)) + rmin;
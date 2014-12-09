function stimulus = InitStim(stimulus,myscreen,categories,dispFig,keepAspectRatio)

%% Image Loading

% make sure widht and height are odd
if iseven(stimulus.widthPix), stimulus.widthPix = stimulus.widthPix-1;end
if iseven(stimulus.heightPix), stimulus.heightPix = stimulus.heightPix-1;end

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
      stimulus.raw{i} = loadNormalizedEqImages(fullfile(stimulus.imageDirMain,categories{i}),'width',stimulus.widthPix,'height',stimulus.heightPix,'dispFig',dispFig,'keepAspectRatio',keepAspectRatio);
      stimulus.p.raw{i} = loadNormalizedEqImages(fullfile(stimulus.imageDirPer,categories{i}),'width',stimulus.p.widthPix,'height',stimulus.p.heightPix,'dispFig',dispFig,'keepAspectRatio',keepAspectRatio);

      % make sure we opened ok
      if isempty(stimulus.raw{i})
        disp(sprintf('(noisecon) Could not load images; %s',categories{i}));
        keyboard
      end
      % get average mag
      stimulus.averageMag = stimulus.averageMag + stimulus.raw{i}.averageMag;
      stimulus.averageDC = stimulus.averageDC +stimulus.raw{i}.averageDC;
      averageN = averageN + 1;
    else
      stimulus.raw{i}.n = 0;
      stimulus.p.raw{i}.n = 0;
    end
  end
  % compute average
  stimulus.averageMag = stimulus.averageMag/averageN;
  stimulus.averageDC = stimulus.averageDC/averageN;
  % and set that we have loaded
  stimulus.imageDir = stimulus.imageDirMain;
  stimulus.p.imageDir = stimulus.imageDirPer;
  stimulus.imagesLoaded = 1;

  disp(sprintf('(noisecon) No file found, generating stimset... \n'));
  
  %% PERIPHERAL IMAGES
  for cat = 1:length(stimulus.p.raw)
      for imgN = 1:stimulus.p.raw{cat}.n
          thisImage = stimulus.p.raw{cat}.halfFourier{imgN};
          %% Get 10-level image
          image = reconstructFromHalfFourier(thisImage);
          
          % Normalize the 10-level image by the equalized histogram
          rmed = .5*255;
          mrmax = 255;
          mrmin = 0;
          % build the normalized PDF
          npdf = normpdf(mrmin:mrmax,rmed,80);
          npdf = npdf / sum(npdf);
          % change the image to match the PDF
          image = (mrmax-mrmin)*histeq(image/255,npdf) + mrmin;
          
          img11 = replaceColors(image,stimulus.colors.nReservedPeripheral);
          stimulus.p.images{cat,1}(:,:,imgN) = img11;
      end
  end
  %% MAIN IMAGES
  for cat = 1:length(stimulus.raw)
      for imgN = 1:stimulus.raw{cat}.n
          thisImage = stimulus.raw{cat}.halfFourier{imgN};
          thisImage.dc = stimulus.averageDC;
          thisImage.mag = stimulus.averageMag;
          
          image = reconstructFromHalfFourier(thisImage);
          
          image = fixBoundaries(image,stimulus.pedestals.maxRange);
          stimulus.tex{cat}(:,:,imgN) = flipud(image);
          
          % Display percent done
          disppercent(calcPercentDone(cat,length(stimulus.raw),imgN,stimulus.raw{cat}.n));
      end
  end
else
  disp(sprintf('(noisecon) Stimulus already initialized'));
end

%% Flip ud or Build Peripheral TExtures

disp(sprintf('(noisecon) Building stimulus array... \n'));
for cat = 1:length(stimulus.p.raw)
    for imgN = 1:stimulus.p.raw{cat}.n                
        i = flipud(stimulus.p.images{cat,1}(:,:,imgN));
        stimulus.p.tex{cat,1}(imgN) = mglCreateTexture(i);
        
        %% Build some 10-level masks (shuffled)
        for m = 2:10
            in = reshape(i(randperm(length(i(:)))),size(i));
            stimulus.p.tex{cat,m}(imgN) = mglCreateTexture(in);
        end
        
        disppercent(calcPercentDone(cat,length(stimulus.p.raw),imgN,stimulus.p.raw{cat}.n));
    end
end

disppercent(inf);

stimulus.trialStart = mglGetSecs;
stimulus.categories = categories;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    imExpand             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rescale an image with any range to be within 0->255 + offset from the
% number of reserved colors. NOTE: this affects the mean luminance!!!

function i3 = imExpand(image)

i3(1,:,:) = image;
i3(2,:,:) = image;
i3(3,:,:) = image;
i3(4,:,:) = 255;
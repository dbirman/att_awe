function stimulus = InitStim(stimulus,myscreen,categories,imageDir,dispFig,keepAspectRatio,saveFile)

if ~isfield(stimulus,'cuenoise'),stimulus.imagesLoaded = 0;end
stimulus.cuenoise = 1;
stimulus.thisTex = [];

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
        disp(sprintf('(noisecon) Could not load images; %s',categories{i}));
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
  disp(sprintf('(noisecon) Stimulus already initialized'));
end

% count how many images
nImages = [];
for i = 1:stimulus.nCategories
  if stimulus.raw{i}.n > 0
    % keep number of images so that we can make scrambles with equal number of images
    nImages(end+1) = stimulus.raw{i}.n;
  end
end

%% Build/Load Image Array

% disp(sprintf('(noisecon) Checking for stimulus file.'));
% if ~exist(saveFile,'file')
%     loaded = 0;
% else
%     loaded = 1;
%     disp(sprintf('(noisecon) File found... loading...'));
%     load(saveFile);
%     disp(sprintf('(noisecon) File found... loaded.'));
% end
loaded = 0;

% note: stimulus.tex{GENDER,RANGE}(:,:,#)
stimulus.LUT = [];

if ~loaded
    disp(sprintf('(noisecon) No file found, generating stimset... '));
    for cat = 1:length(stimulus.raw)
        for imgN = 1:stimulus.raw{cat}.n
            %% Get IMAGE
            thisImage = stimulus.raw{cat}.halfFourier{imgN};
            thisImage.dc = stimulus.averageDC;
            thisImage.mag = stimulus.averageMag;

            %% Get 10-level image
            image = reconstructFromHalfFourier(thisImage);
            img11 = replaceColors(image,stimulus.colors.nReservedPeripheral);            
            stimSave.p.tex{cat,1}(:,:,imgN) = img11;
%             imwrite(img11/255,sprintf('/Users/dan/proj/att_awe/images/output/%s/image.tif',num2str(cat)),'tiff');
            
            % image is currently the image, with the full range (whatever
            % that happens to be)
            
            image = fixBoundaries(image,stimulus.pedestals.maxRange);
            stimSave.tex{cat}(:,:,imgN) = image;
                    
            % Display percent done
            disppercent(calcPercentDone(cat,length(stimulus.raw),imgN,stimulus.raw{cat}.n));
        end
    end
end

%% Save

% if ~loaded
%     disp(sprintf('(noisecon) No file found... saving...'));
%     save(saveFile,'stimSave');
%     disp(sprintf('(noisecon) No file found... saved.'));
% end

%% Build Texture array

disp(sprintf('(noisecon) Building stimulus array... '));
for cat = 1:length(stimulus.raw)
    for imgN = 1:stimulus.raw{cat}.n                
        i = flipud(stimSave.p.tex{cat,1}(:,:,imgN));
        stimulus.p.tex{cat,1}(imgN) = mglCreateTexture(i);
        
        %% Build some 10-level masks (shuffled)
        for m = 2:10
            in = reshape(i(randperm(length(i(:)))),size(i));
            stimulus.p.tex{cat,m}(imgN) = mglCreateTexture(in);
        end
        
        for range = 1:stimulus.pedestals.maxRange
            i = stimSave.tex{cat}(:,:,imgN);
            stimulus.tex{cat}(:,:,imgN) = flipud(i);
        end
        % Display percent done
        disppercent(calcPercentDone(cat,length(stimulus.raw),imgN,stimulus.raw{cat}.n));
    end
end

clear stimSave

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
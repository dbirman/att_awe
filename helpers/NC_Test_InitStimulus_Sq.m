function stimulus = NC_Test_InitStimulus_Sq(stimulus,myscreen,categories,imageDir,dispFig,keepAspectRatio,saveFile)

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

disp(sprintf('(noisecon) Checking for stimulus file.'));
if ~exist(saveFile,'file')
    loaded = 0;
else
    loaded = 1;
    disp(sprintf('(noisecon) File found... loading...'));
    load(saveFile);
    disp(sprintf('(noisecon) File found... loaded.'));
end

% note: stimulus.tex{GENDER,NOISE,NOISE#}(#)
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
            
            stimSave.p.tex{cat,1}(:,:,imgN) = uint8(img11);

            %% Build some 10-level masks (shuffled)
            for m = 2:10
                stimSave.p.tex{cat,m}(:,:,imgN) = uint8(reshape(img11(randperm(length(img11(:)))),size(img11)));
            end

            for exemplar = 1:stimulus.nExemplar
                % Generate a new phase mask
                p_mask = rand(size(stimulus.raw{1}.halfFourier{1}.phase))*2*pi;
                for N = 1:length(stimulus.pedestals.noise)
                    cN = stimulus.pedestals.noise(N);
                    %% NOISE
                    
                    % Use the phase mask to generate a spatially-correlated
                    % noise image
                    exemImg = pinkNoise(image,thisImage,cN,p_mask);
                    exemImg = fixBoundaries(exemImg);
%                     exemImg = replaceColors(exemImg,(256 - stimulus.colors.nReservedColors)/4);
                    % Save Image
                    stimSave.tex{cat,N,exemplar}(:,:,imgN) = uint8(exemImg);
                end
            end
            % Display percent done
            disppercent(calcPercentDone(cat,length(stimulus.raw),imgN,stimulus.raw{cat}.n));
        end
    end
end

%% Save

if ~loaded
    disp(sprintf('(noisecon) No file found... saving...'));
    save(saveFile,'stimSave');
    disp(sprintf('(noisecon) No file found... saved.'));
end

%% Build Texture array

disp(sprintf('(noisecon) Building texture array... '));
for cat = 1:length(stimulus.raw)
    for imgN = 1:stimulus.raw{cat}.n                
        i = imExpand(stimSave.p.tex{cat,1}(:,:,imgN));
        stimulus.p.tex{cat,1}(imgN) = mglCreateTexture(contrastNormalize(flipud(i)));
        
        %% Build some 10-level masks (shuffled)
        for m = 2:10
            i = imExpand(stimSave.p.tex{cat,m}(:,:,imgN));
            stimulus.p.tex{cat,m}(imgN) = mglCreateTexture(contrastNormalize(flipud(i)));
        end
        
        for N = 1:length(stimulus.pedestals.noise)
            for exemplar = 1:stimulus.nExemplar
                % Add the image to the texture set
                i = imExpand(stimSave.tex{cat,N,exemplar}(:,:,imgN));
                stimulus.tex{cat,N,exemplar}(imgN) = mglCreateTexture(contrastNormalize(flipud(i)));
            end
        end
        % Display percent done
        disppercent(calcPercentDone(cat,length(stimulus.raw),imgN,stimulus.raw{cat}.n));
    end
end

clear stimSave

disppercent(inf);

stimulus.trialStart = mglGetSecs;
stimulus.categories = categories;

function image = pinkNoise(img,img_f,K,p_mask)


% Generate white noise
nimg_f = img_f;
nimg_f.phase = p_mask;
noise = reconstructFromHalfFourier(nimg_f);

L0 = mean2(img);

image = L0 + sqrt(K) * (img - L0) + sqrt(1-K) * (noise - L0);

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
%  im = im*0.7;
%  im(im>255) = 255;
%  im(im<0) = 0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %    replaceColors        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
function image2 = replaceColors(image,N)
N = ceil(N);
qts = 1/N:1/N:1;
cutoffs = fliplr(quantile(image(:),qts));
image2 = image;
for i = 1:N
    image2(image<=cutoffs(i)) = N-i;
end

function image = fixBoundaries(image)
global stimulus

image = image * .9;
ma = max(image(:));
image = image - max(ma-255,0);
image(image<=stimulus.colors.nReservedColors) = stimulus.colors.nReservedColors+1;

function i3 = imExpand(image)

% i3 = double(image);
i3(1,:,:) = image;
i3(2,:,:) = image;
i3(3,:,:) = image;
i3(4,:,:) = 255;
% loadNormalizedImages.m
%
%        $Id:$
%      usage: d = loadNormalizedImages(dirname,<width=320>,<height=240>,'dispFig=0')
%         by: justin gardner, daniel birman
%       date: 07/08/11, 11/17/14
%    purpose: loads a set of images in which all the images
%             get the mean amplitude spectrum but their individual
%             phase spectrums
%
function d = loadNormalizedEqImages(dirname,varargin)

if nargin == 0
    help loadNormalizedImages;
    return
end

% get variable arguments
width=[];
height=[];
dispFig=[];
getArgs(varargin,{'height=250','width=250','dispFig=0'});

% check directory
d = [];
if ~isdir(dirname)
    disp(sprintf('(loadNormalizedImages) Could not find directory %s',dirname));
    return
end

% size that image will be resampled to
d.width = width;
d.height = height;

% get a listing of directory
d.dirName = dirname;
d.dir = dir(dirname);
d.n = 0;

mask_size = [100 100];
mask = repmat([zeros(1,75),ones(1,100),zeros(1,75)],100,1);
mask = [zeros(75,250);mask;zeros(75,250)];
% mask = ones(mask_size(1),mask_size(2));

% load each image
if dispFig,smartfig('loadNormalizedImages','reuse');end
disppercent(-inf,sprintf('(loadNormalizedImages) Loading images for %s',dirname));
d.im = zeros(width,height,length(d.dir));
d.averageMag = 0;
d.averageDC = 0;
for i = 1:length(d.dir)
    % get filename
    thisFilename = fullfile(d.dirName,d.dir(i).name);
    % and load if it exists
    if isfile(thisFilename) && (thisFilename(1) ~= '.') && ~isempty(imformats(getext(thisFilename)))
        d.n = d.n + 1;
        % read the image
        [im m alpha] = imread(thisFilename);
        % mask the image
        im = reshape(im(mask==1),mask_size(1),mask_size(2));
%         im = histeq(im);
        % normalize to grayscale and same width height
        im = imageNormalize(im,d.width,d.height,alpha);
        if dispFig,clf;imagesc(im);drawnow;colormap(gray);axis equal; axis off;end
        % save
        d.im(1:width,1:height,d.n) = im;
        d.filenames{d.n} = thisFilename;
        % get its half fourier image
        d.halfFourier{d.n} = getHalfFourier(d.im(:,:,d.n));
        d.averageMag = d.averageMag + d.halfFourier{d.n}.mag;
        d.averageDC = d.averageDC + d.halfFourier{d.n}.dc;
    end
    disppercent(i/length(d.dir));
end
disppercent(inf);
d.im = d.im(:,:,1:d.n);

% now get average magnitude
d.averageMag = d.averageMag/d.n;
d.averageDC = d.averageDC/d.n;


%%%%%%%%%%%%%%%%%%%%%%%%
%    imageNormalize    %
%%%%%%%%%%%%%%%%%%%%%%%%
function im = imageNormalize(im,width,height,alpha)

if ieNotDefined('alpha'),alpha = 255*ones(size(im(:,:,1)));end

% check image dimensions
if ~isequal(size(im(:,:,1)),size(alpha))
    disp(sprintf('(sigdetect:imageNormalize) Alpha image size does not match. Ignoring alpha'));
    alpha = 255*ones(size(im(:,:,1)));
end

% get image dimensions
imdim = size(im);

% first convert to grayscale
if length(imdim > 2)
    im = mean(im,3);
end

% apply alpha (make background gray)
grayvalue = 127;
im = im.*(double(alpha)/255)+grayvalue*(255-double(alpha))/255;

% now resample to the same dimensions
if ~isempty(width) && ~isempty(height)
    [x y] = meshgrid(0:1/(imdim(2)-1):1,0:1/(imdim(1)-1):1);
    [xi yi] = meshgrid(0:1/(height-1):1,0:1/(width-1):1);
    im = interp2(x,y,im,xi,yi,'bilinear');
end

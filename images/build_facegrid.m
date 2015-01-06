function out = build_facegrid()

% This function displays the contrast/noise grid for a random face from the
% specified dataset, using the exact same procedure that cue_noisecon
% employs.

cFolder = '~/proj/att_awe/images/real_faces2_agem/';

gen = {'m/','f/'};

loc = [cFolder gen{2}];

%% Load a random image from loc folder 

imgFiles = dir(strcat(loc,'*jpg'));
oneof = imgFiles(randi(length(imgFiles)));

img = imread([loc oneof.name]);

mask_size = [100 100];
mask = repmat([zeros(1,75),ones(1,100),zeros(1,75)],100,1);
mask = [zeros(90,250);mask;zeros(60,250)];

img = reshape(img(mask==1),mask_size(1),mask_size(2));

%% Setup the contrast and noise increments

% these are from cue_noisecon on 1/6/14
contrast = [ .075 .2 .35 .6 1 ];
noise = [ .119 .269 .5 .731 .881 ];
noise = 1-noise;

%% Build the contrast/noise grid images

out = zeros(size(img,1)-1,size(img,2)-1,length(contrast),length(noise));

mi = 0;
ma = 255;

npdf = normpdf(mi:ma,(ma-mi)/2,50);
i_f = getHalfFourier(img);
p_mask = rand(1,length(i_f.phase))*2*pi;

ctemp = zeros(size(img,1)-1,size(img,2)-1,length(contrast));
% first loop to build the noise images
for i = 1:length(contrast)
    ni = noise(i);
    
    % Add NOISE
    oimg = pinkNoise(i_f,ni,p_mask);
    % rescale to 0-255
    oimg = oimg + abs(min(oimg(:)));
    oimg = oimg * 255 / max(oimg(:));
    
    ctemp(:,:,i) = oimg;
end

% now loop to build the contrast images and add to out
for c = 1:length(contrast)
    ci = contrast(c);
    
    for n = 1:length(noise)
        % get the current contrast image
        cimg = ctemp(:,:,n);

        % scale by the ratio
        cpdf = scalePdf(npdf,mi:ma,ci);
        % adjust contrast
        fimg = histeq(cimg/255,cpdf);
    
        % add to out
        out(:,:,c,n) = fimg;
    end
end

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

function image = pinkNoise(img_f,K,p_mask)

img = reconstructFromHalfFourier(img_f);
% Get the image range
ma = 255;
mi = 0;

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

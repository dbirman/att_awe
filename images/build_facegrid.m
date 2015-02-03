function out = build_facegrid(contrast,noise,p_mask)

% This function displays the contrast/noise grid for a random face from the
% specified dataset, using the exact same procedure that cue_noisecon
% employs.

cFolder = '~/proj/att_awe/images/brazil_faces/';

gen = {'m/','f/'};

loc = [cFolder gen{2}];

%% Load a random image from loc folder 

imgFiles = dir(strcat(loc,'*jpg'));
oneof = imgFiles(31);

img = imread([loc oneof.name]);

% mask_size = [100 100];
% mask = repmat([zeros(1,75),ones(1,100),zeros(1,75)],100,1);
% mask = [zeros(90,250);mask;zeros(60,250)];
% 
% img = reshape(img(mask==1),mask_size(1),mask_size(2));

%% Setup the contrast and noise increments

% these are from cue_noisecon on 1/6/14
% noise = 1-noise;

%% Build the contrast/noise grid images

out = zeros(size(img,1),size(img,2)-1,length(contrast),length(noise));

mi = 0;
ma = 255;

npdf = normpdf(mi:ma,(ma-mi)/2,55);
i_f = getHalfFourier(img);

ctemp = zeros(size(img,1),size(img,2)-1,length(contrast));
% first loop to build the noise images
for i = 1:length(noise)
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


%% test func

figure
hold on
plot(0:255,npdf)
for x = fliplr(.01:.3:.91)
s = scalePdf(npdf,0:255,x);
plot(0:255/(length(s)-1):255,s)
end
axis([0 255 .0005 .008])
legend({'Full Contrast','90%','60%','30%','1%'})
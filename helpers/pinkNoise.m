
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

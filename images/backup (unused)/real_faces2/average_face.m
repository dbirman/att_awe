%% Load files

cFolder = '~/proj/att_awe/images/real_faces2_agem/';

gen = {'m/','f/'};

%% Load face data
data.fimages = zeros(100,100,95);
data.mimages = zeros(100,100,94);

mask_size = [100 100];
mask = repmat([zeros(1,75),ones(1,100),zeros(1,75)],100,1);
mask = [zeros(90,250);mask;zeros(60,250)];

gender = gen{1};
loc = [cFolder gender];
imgFiles = dir(strcat(loc,'*jpg'));
for imgNum = 1:length(imgFiles)
    im = imread([loc imgFiles(imgNum).name]);
    im = reshape(im(mask==1),mask_size(1),mask_size(2));
    data.mimages(:,:,imgNum) = im;
end

gender = gen{2};
loc = [cFolder gender];
imgFiles = dir(strcat(loc,'*jpg'));
for imgNum = 1:length(imgFiles)
    im = imread([loc imgFiles(imgNum).name]);
    im = reshape(im(mask==1),mask_size(1),mask_size(2));
    data.fimages(:,:,imgNum) = im;
end

%% Plot
figure
subplot(121)
colormap('gray');
imagesc(mean(data.mimages,3));
subplot(122)
colormap('gray');
imagesc(mean(data.fimages,3));
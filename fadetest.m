function fadetest()
%% Init a test screen

mglOpen
mglVisualAngleCoordinates(57,[13 12]);
mglSetParam('visualAngleSquarePixels',0,1);
mglClearScreen(.5)
mglFlush

%% Load some images

files.raw = loadNormEqImages2('/Users/dan/proj/att_awe/images/all_faces2/f/','width',399,'height',399,'dispFig',0,'keepAspectratio',1);

%% Set PArams

npdf = normpdf(0:255,125.5,50);
npdf = scalePdf(npdf,0:255,1);
% scale by the ratio
% change the image to match the PDF
gaussianWin = mglMakeSmoothBorderMask(11,12,8.5,2);
win = 255-255*(gaussianWin);
mask = ones(size(win,1),size(win,2),4)*.5*255;
mask(:,:,4) = win;
maskTex = mglCreateTexture(mask);
%% Put an image on the screen
img = flipud(files.raw.im(:,:,1));
img2 = (255-0)*histeq(img/255,npdf) + 0;

imgTex = mglCreateTexture(img2);
mglBltTexture(imgTex,[0 0 4 4]);
mglBltTexture(maskTex,[0 0 5 5]);
mglFlush

%% Fade into a low contrast version

steps = 100;
goal = .1;

imgs = {};
for i = 1:steps
    img2 = (255-0)*histeq(img/255,scalePdf(npdf,0:255,(goal-i*((1-goal)/steps))+goal)) + 0;
    imgs{i} = mglCreateTexture(img2);
end

%% Go!
% mglBltTexture(imgTex,[0 0 4 4]);
% mglBltTexture(maskTex,[0 0 5 5]);
% mglFlush
for i = 1:steps
    mglWaitSecs(1);
    mglBltTexture(imgs{i},[0 0 4 4]);
    mglBltTexture(maskTex,[0 0 5 5]);
    mglFlush
%     disp('Going');
end

%% test

%  mglOpen;
%   mglVisualAngleCoordinates(57,[52 32.7]);
%   g = mglMakeGrating(10,10,1.5,90,0);
%   m = mglMakeSmoothBorderMask(10,10,8,1.5,3);
%   g = g.*m;
%   g = 255*(g+1)/2;
%   tex = mglCreateTexture(g);
%   mglClearScreen(0.5);
%   mglBltTexture(tex,[0 0]);
%   mglFlush;
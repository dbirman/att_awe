%%
basepdf = normpdf(0:255,125.5,75);

image = imread('/Users/dan/proj/att_awe/images/all_faces2/f/face3.bmp');
image = mean(image,3);

image2 = imread('/Users/dan/proj/att_awe/images/all_faces2/f/face2.bmp');
image2 = mean(image2,3);

%%
npdf = scalePdf(basepdf,0:255,.75);

img = 255*histeq(image/255,npdf);
img2 = 255*histeq(image2/255,npdf);

figure(1)
colormap('gray')

%%

for i = 0:100 % i ist he percent of image2 we use at each pixel
    imagesc(img*(100-i)/100 + img2*i/100,[0 255]);
    pause(.01)
end

for i = fliplr(0:100) % i ist he percent of image2 we use at each pixel
    imagesc(img*(100-i)/100 + img2*i/100,[0 255]);
    pause(.01)
end
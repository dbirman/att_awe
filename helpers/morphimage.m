%morphimage morphs an image (grayscale) using a morphing matrix.
%	MorphedImage = morphimage(Image, morphfield) morphs Image.
%		
%	MorphedImage --> Result of the morphing process
%
%   Image        --> Original grayscale image
%   morphfield   --> morphing matrix of size m,n,2 where
%                       :,:,1 --> morphing matrix for x-coordiantes
%                       :,:,2 --> morphing matrix for y-coordiantes
%
%example for a grayscale image: 
%Read a grayscale image
%   image = imread('/Users/dan/proj/att_awe/images/all_faces2/f/face1.bmp'); 
%   image = mean(image,3);
%   imagesize = size(image);
% Creating a random morphield
%   morphfield = 8*imresize(rand(5,5,2)-0.5,imagesize,'bilinear');
% % Morph the original image
%   MorphedImage = morphimage(image,morphfield);
% % Show original image
	figure(1); imagesc(image,[0 255]); colormap gray;
% % Show morphed image
	figure(2); imagesc(MorphedImage,[0 255]); colormap gray;
% % Show morphing matrix for x-coordinates
	figure(3); imagesc(morphfield(:,:,1)); colorbar;
% % Show morphing matrix for y-coordinates
	figure(4); imagesc(morphfield(:,:,2)); colorbar;


function res = morphimage(orig,morphfield)

%get size of the image
origsize = size(orig);
%create an empty image
res = zeros(origsize);

%compute the new grayvalues for every pixel
for cnt1 = 1:1:origsize(1)
    for cnt2 = 1:1:origsize(2)
        %get the morphing-vector for the current pixel
        xn = cnt2-morphfield(cnt1,cnt2,1);
        yn = cnt1-morphfield(cnt1,cnt2,2);
        if xn < 1
            xn = 1;
        end
        if xn > origsize(2)
            xn = origsize(2);
        end
        if yn < 1
            yn = 1;
        end
        if yn > origsize(1)
            yn = origsize(1);
        end
        %compute the pixel coordinates (x1,x2,y1,y2) and weighting 
        %(fx1,fx2,fy1,fy2) which are used to compute the new 
        %pixel grayvalue
        x1 = floor(xn);
        x2 = ceil(xn);
        fx1 = x2-xn;
        fx2 = 1-fx1;
        y1 = floor(yn);
        y2 = ceil(yn);
        fy1 = y2-yn;
        fy2 = 1-fy1;
        %compute the new grayvalue of the current pixel
        res(cnt1,cnt2) = fx1*fy1*orig(y1,x1)+fy1*fx2*orig(y1,x2)+fy2*fx1*orig(y2,x1)+fy2*fx2*orig(y2,x2);
    end
end
    

        
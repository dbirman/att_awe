%% get phase random image

random = rand(size(img));

i_rnd = fft2(random);
i_fft = fft2(img);

pimg = ifft2(abs(i_fft).*exp(sqrt(-1)*angle(i_rnd)));
pimg = real(pimg);
pimg = (pimg + min(pimg(:)));
pimg = pimg * 255 / max(pimg(:));
%% generate % contrasts and phase randomization
cons = 0:.1:1;
phs = 0:.1:1;
figure; hold on
w = size(img,1); h = size(img,2);
for ci = 1:11
    con = cons(ci);
    for pi = 1:11
        ph = phs(pi);
        
        % create new image as phase blend
        cimg = ph*pimg + (1-ph)*double(img);
        % adjust contrast
        cimg = cimg-127.5;
        cimg = cimg * con;
        cimg = cimg+127.5;
        
%         oimg(ci,pi,:,:) = cimg;
        
%         subplot(11,11,(ci-1)*11+pi);
        xidx = (ci-1)*w+1 : ci*w;
        yidx = (pi-1)*h+1 : pi*h;
        imagesc(xidx,yidx,flipud(cimg)); colormap('gray');
%         disp(caxis);
%         caxis([0 255]);
%         set(gca,'Visible','off');
    end
end

axis square
set(gca,'XTick',[],'YTick',[]);
drawPublishAxis;
close all hidden
p_mask = rand(1,15536)*2*pi;

creps = {.5,0:.01:1};
nreps = {0:.01:1,.5};

p_mask = rand(1,15536)*2*pi;
% Load an image and its grid
out = build_facegrid(.1:.1:.9,.1:.1:.9,p_mask);

f = figure

% Build all the plots
for i = 1:size(out,3)
    for j = 1:size(out,4)
        subplot(size(out,3),size(out,4),(i-1)*5+j);
        imagesc(out(:,:,i,j),[0 1]);
        colormap('gray');
        axis([1 100 1 100]);
    end
end

% get rid of the axes
set(findobj(gcf, 'type','axes'), 'Visible','off')

% print
print(f,'-dpdf','~/proj/att_awe/images/grid_image.pdf');
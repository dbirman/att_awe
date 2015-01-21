close all hidden

contrast = [ .075 .2 .35 .6 1 ];
% lognoise -2:2
noise = [ 0.119202922022118         0.268941421369995                       0.5         0.731058578630005         0.880797077977882 ];
snr = noise ./ (1-noise);
% Load an image and its grid
out = build_facegrid(contrast,noise);

% OUT(:,:,CONTRAST,NOISE)

f = figure
colormap('gray');
axis([1 100 1 100]);

% Build all the plots
for i = 1:size(out,3)
    for j = 1:size(out,4)
        imagesc(out(:,:,i,j),[0 1]);
        set(findobj(gcf, 'type','axes'), 'Visible','off')
        print(f,'-dpdf',fullfile('~/proj/att_awe/images/built/',sprintf('con%snoi%s_grid_image.pdf',num2str(i),num2str(j))));

    end
end

close all

% get rid of the axes

% print
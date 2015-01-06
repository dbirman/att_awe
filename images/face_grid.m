close all hidden

%out = build_facegrid;

f = figure

for i = 1:size(out,3)
    for j = 1:size(out,4)
        subplot(size(out,3),size(out,4),(i-1)*5+j);
        imagesc(out(:,:,i,j),[0 1]);
        colormap('gray');
        axis([1 100 1 100]);
    end
end

set(findobj(gcf, 'type','axes'), 'Visible','off')

print(f,'-dpdf','~/proj/att_awe/images/grid_image.pdf');
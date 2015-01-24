
creps = {.5,0:.01:1};
nreps = {0:.01:1,.5};

% p_mask = rand(1,15536)*2*pi;
for run = 1:length(creps)
    contrast = creps{run};
    noise = nreps{run};
    snr = noise ./ (1-noise);
    % Load an image and its grid
    out = build_facegrid(contrast,noise,p_mask);
    
    % OUT(:,:,CONTRAST,NOISE)
    
    % f = figure('visible','off');
    % colormap('gray');
    % axis([1 100 1 100]);
    
    % Build all the plots
    for i = 1:size(out,3)
        for j = 1:size(out,4)
            %         imagesc(out(:,:,i,j),[0 1]);
            %         set(findobj(gcf, 'type','axes'), 'Visible','off')
            fname = fullfile('~/proj/att_awe/images/built/',sprintf('con%03.fnoi%03.f_img.png',contrast(i)*100,100*noise(j)));
            imwrite(out(:,:,i,j),fname,'PNG');
            %         print(f,strcat('-dpdf ',));
            %         print(f,'-dtiff', fname);
        end
    end

end

% get rid of the axes

% print
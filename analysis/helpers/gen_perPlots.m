function gen_perPlots(data,norm)

figure
hold on
means = []; % single, dual noise, dual contrast

normval = 0;

for j = 1:3
    dat = data{j};
    select = dat(logical([dat<.25].*[dat>0]));
    if j == 1
        normval = mean(select);
    end
    if norm
        select = select / normval;
    end
    means(j) = mean(select);
    stds(j) = std(select)/sqrt(length(select));
end


barwitherr(stds,means);
set(gca,'XTickLabel',{'Single','Dual: Noise','Dual: Contrast'})
set(gca,'XTick',[1 2 3]);
title('Gender');
ylabel('Stimulus Onset Asynchrony (ms)');
barmap=[0.7 0.7 0.7; 0.05 .45 0.1]; %[0.7 0.7 0.7] is grey, [ 0.05 .45 0.1] is green
colormap(barmap)

normed = {'','NORM'};
try
    if isempty(mglGetSID)
        disp('Set SID before printing!');
    end
    print(gcf,'-dpdf',fullfile('~/proj/att_awe/analysis/',mglGetSID,sprintf('figures/GenderDiscriminationFunction%s',normed{norm+1})));
catch
    warning('Print failed...');
end
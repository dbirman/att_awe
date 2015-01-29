function gen_discPlots(plotting,plottingsd,stimulus)


dispTexts = {'Noise','Contrast'};
colorOpts = [1 0 0
             0 1 0];
typePs = {'SnR','contrast'};
stimulus.pedestals.lognoise = log(stimulus.pedestals.noise./(1-stimulus.pedestals.noise));
stimulus.pedestals.SnR = stimulus.pedestals.noise./(1-stimulus.pedestals.noise);

for num = 1:2
    dispText = dispTexts{num};
    color = colorOpts(num,:);
    typeP = typePs{num};
    % Discrimination function plots
    figure
    hold on
    title(sprintf('%s',dispText));
% % %     Plotting(CUES 1/4, PED 1:3, DUAL, RUNTYPE noise/con);
% % % % %     if num == 1
% % % % %         dat = plotting(:,:,:,1);
% % % % %         dat = dat + repmat(stimulus.pedestals.noise(2:4),2,1,2);
% % % % %         dat = dat ./ (1-dat);
% % % % %         dat = dat ./ repmat(stimulus.pedestals.SnR(2:4),2,1,2);
% % % % %         plotting(:,:,:,1) = dat;
% % % % %     end

    errorbar(stimulus.pedestals.(typeP)(2:4),plotting(1,:,1,num),plottingsd(1,:,1,num),'-','Color',color); % FOCAL, PEDS 1:3, SINGLE
    errorbar(stimulus.pedestals.(typeP)(2:4),plotting(2,:,1,num),plottingsd(2,:,1,num),'--','Color',color); % DIST, SINGLE
    errorbar(stimulus.pedestals.(typeP)(2:4),plotting(1,:,2,num),plottingsd(1,:,2,num),'-','Color',.5*color); % FOCAL, DUAL
    errorbar(stimulus.pedestals.(typeP)(2:4),plotting(2,:,2,num),plottingsd(2,:,2,num),'--','Color',.5*color); % DIST, DUAL
    
    if num == 1
        a = axis;
        axis([stimulus.pedestals.(typeP)(2)-.05 stimulus.pedestals.(typeP)(4)+.05 0 a(4)]);
        xlabel('Noise Pedestals (SnR units)');
        ylabel('Delta Noise at Threshold (arbitrary units)');
    else
        a = axis;
        axis([stimulus.pedestals.(typeP)(2)-.025 stimulus.pedestals.(typeP)(4)+.025 0 a(4)]);
        xlabel('Contrast Pedestals (%)');
        ylabel('Delta Contrast at Threshold (%)');
%         set(gca,'xscale','log');
    end
    legend('Focal, Single Task','Distributed, Single Task','Focal, Dual Task','Distributed, Dual Task');
    normed = {'', 'NORM'};
    try
        if isempty(mglGetSID)
            disp('Set SID before printing!');
        end
        print(gcf,'-dpdf',fullfile('~/proj/att_awe/analysis/',mglGetSID,sprintf('figures/%sDiscriminationFunction%s',dispText,normed{norm+1})));
    catch
        warning('Print failed...');
    end
    hold off
end

function gen_discFuncs(stimulus)

disp('Computing Weibull functions. CAUTION: Check all Weibull functions for accuracy (use check=1 flag)');
check = 0;

dispTexts = {'Noise','Contrast'};
colorOpts = [1 0 0
             0 1 0];
typePs = {'noise','contrast'};
plotting = zeros(2,3,2,2);
for num = 1:2
    dispText = dispTexts{num};
    color = colorOpts(num,:);
    typeP = typePs{num};
    
    for cues = 1:2
        for ped = 1:3
            try
                if check
                    out1 = doStaircase('threshold',stimulus.staircase{num,cues,ped},'dispFig',1,'type','weibull'); % noise, 1 cue, lowest
                    keyboard
                else
                    out1 = doStaircase('threshold',stimulus.staircase{num,cues,ped},'type','weibull'); % noise, 1 cue, lowest
                end
                plotting(cues,ped,1,num) = out1.threshold;
            catch
                plotting(cues,ped,1,num) = -1;
            end
            try
                if check
                    out2 = doStaircase('threshold',stimulus.dualstaircase{num,cues,ped},'dispFig',1,'type','weibull'); % noise, 1 cue, lowest
                    keyboard
                else
                    out2 = doStaircase('threshold',stimulus.dualstaircase{num,cues,ped},'type','weibull'); % noise, 1 cue, lowest
                end
                plotting(cues,ped,2,num) = out2.threshold;
            catch
                plotting(cues,ped,2,num) = -1;
            end
        end
    end
    % Discrimination function plots
    figure
    hold on
    title(sprintf('%s',dispText));
    plot(stimulus.pedestals.(typeP)(2:4),plotting(1,:,1,num),'-','Color',color);
    plot(stimulus.pedestals.(typeP)(2:4),plotting(2,:,1,num),'--','Color',color);
    plot(stimulus.pedestals.(typeP)(2:4),plotting(1,:,2,num),'-','Color',.5*color);
    plot(stimulus.pedestals.(typeP)(2:4),plotting(2,:,2,num),'--','Color',.5*color);
    if num == 1
        axis([stimulus.pedestals.(typeP)(2)-.1 stimulus.pedestals.(typeP)(4)+.1 0 .8]);
    else
        axis([stimulus.pedestals.(typeP)(2)-.05 stimulus.pedestals.(typeP)(4)+.05 0 .6]);
    end
    legend('Focal, Single Task','Distributed, Single Task','Focal, Dual Task','Distributed, Dual Task');
    try
        print(gcf,'-dpdf',sprintf('~/proj/att_awe/analysis/figures/%sDiscriminationFunction',dispText));
    catch
        warning('Print failed...');
    end
    hold off
end

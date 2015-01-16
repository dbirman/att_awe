function plotting = gen_discFuncs(stimulus)

disp('Computing Weibull functions. CAUTION: Check all Weibull functions for accuracy (use check=1 flag)');
check = 0;

dispTexts = {'Noise','Contrast'};
colorOpts = [1 0 0
             0 1 0];
typePs = {'lognoise','contrast'};
stimulus.pedestals.lognoise = log(stimulus.pedestals.noise./(1-stimulus.pedestals.noise));
plotting = zeros(2,3,2,2);
plottingsd = zeros(2,3,2,2);
for num = 1:2
    dispText = dispTexts{num};
    color = colorOpts(num,:);
    typeP = typePs{num};
    
    for cues = 1:2
        for ped = 1:3
            try
                if check
%                     out1 = doStaircase('threshold',stimulus.staircase{num,cues,ped},'dispFig',1,'type','weibull'); % noise, 1 cue, lowest
%                     keyboard
                else
                    for i = 1:length(stimulus.staircase{num,cues,ped})
                        out{i} = doStaircase('threshold',stimulus.staircase{num,cues,ped}(i),'type','weibull');% noise, 1 cue, lowest
                    end 
                end
                thresh = [];
                for i = 1:length(out)
                    thresh(i) = out{i}.threshold;
                end
                plotting(cues,ped,1,num) = mean(thresh);
                plottingsd(cues,ped,1,num) = std(thresh)/sqrt(length(thresh));
            catch
                plotting(cues,ped,1,num) = -1;
            end
            try
                if check
                    out2 = doStaircase('threshold',stimulus.dualstaircase{num,cues,ped},'dispFig',1,'type','weibull'); % noise, 1 cue, lowest
                    keyboard
                else
                    for i = 1:length(stimulus.dualstaircase{num,cues,ped})
                        out2{i} = doStaircase('threshold',stimulus.dualstaircase{num,cues,ped}(i),'type','weibull');% noise, 1 cue, lowest
                    end
                end
                for i = 1:length(out)
                    thresh(i) = out2{i}.threshold;
                end
                plotting(cues,ped,2,num) = mean(thresh);
                plottingsd(cues,ped,2,num) = std(thresh)/sqrt(length(thresh));
            catch
                plotting(cues,ped,2,num) = -1;
            end
        end
    end
    % Discrimination function plots
    figure
    hold on
    title(sprintf('%s',dispText));
    if num == 1
        plotting = log(plotting);
    end
    errorbar(stimulus.pedestals.(typeP)(2:4),plotting(1,:,1,num),plottingsd(1,:,1,num),'-','Color',color); % FOCAL, PEDS 1:3, SINGLE
    errorbar(stimulus.pedestals.(typeP)(2:4),plotting(2,:,1,num),plottingsd(2,:,1,num),'--','Color',color); % DIST, SINGLE
    errorbar(stimulus.pedestals.(typeP)(2:4),plotting(1,:,2,num),plottingsd(1,:,2,num),'-','Color',.5*color); % FOCAL, DUAL
    errorbar(stimulus.pedestals.(typeP)(2:4),plotting(2,:,2,num),plottingsd(2,:,2,num),'--','Color',.5*color); % DIST, DUAL
    end
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

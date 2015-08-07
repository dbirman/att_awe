function [plotting, fits, allfit] = cohCon_nocatPerf(stimulus,check)

plotting = cell(2,1);
fits = cell(2,1);
allfit = cell(2,1);
try
    for task = 1:2
        for ped = 1
            if check
            end
            try
                for i = 1:length(stimulus.nocatchs.staircase{task,ped})
                    if stimulus.nocatchs.staircase{task,ped}(i).trialNum>0
                        out = doStaircase('threshold',stimulus.nocatchs.staircase{task,ped}(i),'type','weibull'); % noise, 1 cue, lowest
                        plotting{task,ped} = [plotting{task,ped} out.threshold];
                        fits{task,ped} = [fits{task,ped}; out.fit.fitparams];
                    end
                end
            catch
                plotting{task,ped} = -1;
            end
        end
        allfit{task} = doStaircase('threshold',stimulus.nocatchs.staircase{task,ped},'type','weibull');
    end
catch
    disp('(cohCon_catPerf) Failed.');
end
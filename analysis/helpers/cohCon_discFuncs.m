function [plotting, fits, allfit] = cohCon_discFuncs(stimulus,check)

disp('Computing Weibull functions. CAUTION: Check all Weibull functions for accuracy (use check=1 flag)');

plotting = cell(2,1);
fits = cell(2,1);
allfit = cell(2,1);
for task = 1:2
    
    for ped = 1
        if check
            for i = 1:length(stimulus.staircase{task,ped})
                if stimulus.staircase{task,ped}(i).trialNum>0
                    out = doStaircase('threshold',stimulus.staircase{task,ped}(i),'type','weibull','dispFig=1'); % noise, 1 cue, lowest
                    plotting{task,ped} = [plotting{task,ped} out.threshold];
                    fits{task,ped} = [fits{task,ped}; out.fit.fitparams];
                    keyboard
                end
            end
        else
            for i = 1:length(stimulus.staircase{task,ped})
                if stimulus.staircase{task,ped}(i).trialNum>0
                    out = doStaircase('threshold',stimulus.staircase{task,ped}(i),'type','weibull'); % noise, 1 cue, lowest
                    plotting{task,ped} = [plotting{task,ped} out.threshold];
                    fits{task,ped} = [fits{task,ped}; out.fit.fitparams];
                end
            end
        end
        allfit{task} = doStaircase('threshold',stimulus.staircase{task,ped},'type','weibull');
    end
end

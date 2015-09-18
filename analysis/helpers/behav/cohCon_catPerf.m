function [plotting, fits, allfit] = cohCon_catPerf(stimulus,check)

plotting = cell(2,1);
fits = cell(2,1);
allfit = cell(2,1);
for task = 1:2
    for ped = 1
        if check
        end
        for i = 1:length(stimulus.stairCatch{task,ped})
            if stimulus.stairCatch{task,ped}(i).trialNum>0
                out = doStaircase('threshold',stimulus.stairCatch{task,ped}(i),'type','weibull'); % noise, 1 cue, lowest
                if out.threshold > 1
                    out = doStaircase('threshold',stimulus.stairCatch{task,ped}(i));
                    t = out.meanOfAllReversals;
                else
                    t = out.threshold;
                end
                plotting{task,ped} = [plotting{task,ped} t];
                if isfield(out,'fit')
                    fits{task,ped} = [fits{task,ped}; out.fit.fitparams];
                end
            end
        end
        allfit{task} = doStaircase('threshold',stimulus.stairCatch{task,ped},'type','weibull');
    end
end
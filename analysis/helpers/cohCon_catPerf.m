function plotting = cohCon_catPerf(stimulus,check)

plotting = cell(2,1);
try
    for task = 1:2
        for ped = 1
            if check
            end
            try
                for i = 1:length(stimulus.stairCatch{task,ped})
                    if stimulus.stairCatch{task,ped}(i).trialNum>0
                        out = doStaircase('threshold',stimulus.stairCatch{task,ped}(i),'type','weibull'); % noise, 1 cue, lowest
                        plotting{task,ped} = [plotting{task,ped} out.threshold];
                    end
                end
            catch
                plotting{task,ped} = -1;
            end
        end
    end
catch
    disp('(cohCon_catPerf) Failed.');
end
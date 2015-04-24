function plotting = cohCon_discFuncs(stimulus,check)

disp('Computing Weibull functions. CAUTION: Check all Weibull functions for accuracy (use check=1 flag)');

plotting = cell(2,1);
try

    for task = 1:2

        for ped = 1
            if check
                try
                    for i = 1:length(stimulus.staircase{task,ped})
                        if stimulus.staircase{task,ped}(i).trialNum>0
                            out = doStaircase('threshold',stimulus.staircase{task,ped}(i),'type','weibull','dispFig=1'); % noise, 1 cue, lowest
                            plotting{task,ped} = [plotting{task,ped} out.threshold];
                            keyboard
                        end
                    end
                catch
                    plotting{task,ped} = -1;
                end
            else
                try
                    for i = 1:length(stimulus.staircase{task,ped})
                        if stimulus.staircase{task,ped}(i).trialNum>0
                            out = doStaircase('threshold',stimulus.staircase{task,ped}(i),'type','weibull'); % noise, 1 cue, lowest
                            plotting{task,ped} = [plotting{task,ped} out.threshold];
                        end
                    end
                catch
                    plotting{task,ped} = -1;
                end
            end
        end
    end
catch
    disp('(cohCon_discFuncs) Failed.');
end
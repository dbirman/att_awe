function plotting = cohCon_discFuncs(stimulus,check)

disp('Computing Weibull functions. CAUTION: Check all Weibull functions for accuracy (use check=1 flag)');

plotting = cell(2,4);
try

    for task = 1:2

        for ped = 1:4
            if check
            end
            try
                for i = 1:length(stimulus.staircase{task,ped})
                    out = doStaircase('threshold',stimulus.staircase{task,ped}(i),'type','weibull'); % noise, 1 cue, lowest
                    plotting{task,ped} = [plotting{task,ped} out.threshold];
                end
            catch
                plotting{task,ped} = -1;
            end
        end
    end
catch
    disp('(cohCon_discFuncs) Failed.');
end
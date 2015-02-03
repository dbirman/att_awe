function plotting = gen_discFuncs(stimulus,check)

disp('Computing Weibull functions. CAUTION: Check all Weibull functions for accuracy (use check=1 flag)');

plotting = cell(2,3);
    for cues = 1:2
        for ped = 1:3
            try
                if check
%                     out1 = doStaircase('threshold',stimulus.staircase{num,cues,ped},'dispFig',1,'type','weibull'); % noise, 1 cue, lowest
%                     keyboard
                else
                    for i = 1:length(stimulus.staircase{cues,ped})
                        if stimulus.staircase{cues,ped}(i).trialNum > 0
                            out{i} = doStaircase('threshold',stimulus.staircase{cues,ped}(i),'type','weibull');% noise, 1 cue, lowest
                        end
                    end 
                end
                thresh = [];
                for i = 1:length(out)
                    thresh(i) = out{i}.threshold;
                end
                plotting{cues,ped} = thresh;
            catch
                plotting{cues,ped} = -1;
            end
        end
    end

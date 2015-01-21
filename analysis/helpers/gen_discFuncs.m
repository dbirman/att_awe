function plotting = gen_discFuncs(stimulus)

disp('Computing Weibull functions. CAUTION: Check all Weibull functions for accuracy (use check=1 flag)');
check = 0;

plotting = cell(2,3,2,2);
for num = 1:2
    
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
                plotting{cues,ped,1,num} = thresh;
            catch
                plotting{cues,ped,1,num} = -1;
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
                for i = 1:length(out2)
                    thresh2(i) = out2{i}.threshold;
                end
                plotting{cues,ped,2,num} = thresh2;
            catch
                plotting{cues,ped,2,num} = -1;
            end
        end
    end
end
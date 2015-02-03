function staircaseplots( stimulus )
% Depending on the type of the loaded run, generate the staircase figures
stairtype = 'dualstaircase';

dispTexts = {'Noise','Contrast'};
for num = 1:2
    plotting = zeros(2,3);
    
    dispText = dispTexts{num};
    
    try
        figure % this is the 'staircase' figure
        title(sprintf('%s, DUAL Staircase plot (R->G->B high)',dispText));
        hold on
        drawing = {'-r' '-g' '-b'
            '--r' '--g' '--b'};
        for cues = 1:2
            for ped = 1:3
                testV = [];
                % for each cue/ped combo, this pulls out and adds to the test
                % value list. This will be used to display the 'staircases'.
                for i = 1:length(stimulus.(stairtype){num,cues,ped})
                    testV = [testV stimulus.(stairtype){num,cues,ped}(i).testValues];
                end
                plot(testV,drawing{cues,ped});
            end
        end
    catch
        warning('Dual task staircases failed');
    end
    stairtype = 'staircase';
    try
        figure % this is the 'staircase' figure
        title(sprintf('%s, Staircase plot (R->G->B high)',dispText));
        hold on
        drawing = {'-r' '-g' '-b'
            '--r' '--g' '--b'};
        for cues = 1:2
            for ped = 1:3
                testV = [];
                % for each cue/ped combo, this pulls out and adds to the test
                % value list. This will be used to display the 'staircases'.
                for i = 1:length(stimulus.(stairtype){num,cues,ped})
                    testV = [testV stimulus.(stairtype){num,cues,ped}(i).testValues];
                end
                plot(testV,drawing{cues,ped});
            end
        end
        if stimulus.p.staircase(end).trialNum == 0
            % Ignore the last staircase, it was just reset and has no trials
            % (and causes an error)
            usestair = stimulus.p.staircase(1:end-1);
        else
            usestair = stimulus.p.staircase;
        end
    catch
        warning('Single task staircases failed');
    end
end
% Peripheral task
try
    disp('Gender performance during NOISE task');
    doStaircase('threshold',stimulus.p.dualstaircase{1},'dispFig',1,'type','weibull');
    disp('Gender performance during CONTRAST task');
    doStaircase('threshold',stimulus.p.dualstaircase{2},'dispFig',1,'type','weibull');
    disp('Gender performance during SOLO task');
    doStaircase('threshold',usestair,'dispFig',1,'type','weibull');
catch
    warning('Peripheral task staircases failed');
end

end


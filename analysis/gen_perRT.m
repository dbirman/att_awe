function gender = gen_perRT( expHolder )

for i = 1:length(expHolder)
    cExp = expHolder{i};
%     main = cExp{1};
    per = cExp{2};
    rVars = cExp{3}.runVars;
    
    % Peripheral Task Differences
    for t = 2:length(per.trials)
        if per.randVars.respond(t)==0
            type = 3;
        else
            type = per.randVars.tGen(t);
        end
        correct = per.correct(t)+1;
        if isnan(correct), correct = 3; end
        try
            if per.randVars.mainTrialNum(t+1) == per.randVars.mainTrialNum(t)
                order = 2;
            else
                order = 1;
            end
        catch
            % TODO: remove this when all runs include the correct
            % mainTrialNum flag...
            order = randi(2);
        end
        gender(rVars.dual+1,end+1,order,correct,type) = per.reactionTime(t);
        disppercent(calcPercentDone(t,length(per.trials)),'Calculating: ');
    end
    
    % Main task differences
end

end


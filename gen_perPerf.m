function data = gen_perPerf(stimulus)


data = {};
for i = 1:length(stimulus.p.staircase)
    if ~stimulus.p.staircase(i).trialNum == 0
        out = doStaircase('threshold',stimulus.p.staircase(i),'type','weibull');
        data{1}(i) = out.threshold;
    end
end

for dual = 1:2
    cur = stimulus.p.dualstaircase{dual};
    for k = 1:length(cur)
        if ~cur(k).trialNum == 0
            out = doStaircase('threshold',cur(k),'type','weibull');
            data{dual+1}(k) = out.threshold;
        end
    end
end

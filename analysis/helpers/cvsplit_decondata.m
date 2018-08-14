function dd = cvsplit_decondata( dd)

ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
datas = {'cc','time'};

cfold = 1;

for cidx = 1:length(dd.V1.cc.conidxs)
    tidx = [];
    fold = struct;
    cidxtrain = setdiff(1:length(dd.V1.cc.conidxs),cidx);
%     tidxtrain = 1:length(dd.V1.time.conidxs);

    for ri = 1:8
        fold.(ROIs{ri}).train = dd.(ROIs{ri});
        fold.(ROIs{ri}).test = dd.(ROIs{ri});
        
        % remove all but the test from cc
        fs = fields(fold.(ROIs{ri}).test.time);
        for fi = 1:length(fs)
            eval(sprintf('fold.(ROIs{ri}).test.time.%s = fold.(ROIs{ri}).test.time.%s(tidx);',fs{fi},fs{fi}));
        end
        fs = fields(fold.(ROIs{ri}).test.cc);
        for fi = 1:length(fs)
            eval(sprintf('fold.(ROIs{ri}).test.cc.%s = fold.(ROIs{ri}).test.cc.%s(cidx);',fs{fi},fs{fi}));
        end
        for fi = 1:length(fs)
            eval(sprintf('fold.(ROIs{ri}).train.cc.%s = fold.(ROIs{ri}).train.cc.%s(cidxtrain);',fs{fi},fs{fi}));
        end
    end
    dd.(sprintf('fold%i',cfold)) = fold; cfold=cfold+1;
end

for tidx = 1:length(dd.V1.time.conidxs)
    cidx = [];
    fold = struct;
%     cidxtrain = setdiff(1:length(dd.V1.cc.conidxs),cidx);
    tidxtrain = setdiff(1:length(dd.V1.time.conidxs),tidx);

    for ri = 1:8
        fold.(ROIs{ri}).train = dd.(ROIs{ri});
        fold.(ROIs{ri}).test = dd.(ROIs{ri});
        
        % remove all but the test from cc
        fs = fields(fold.(ROIs{ri}).test.cc);
        for fi = 1:length(fs)
            eval(sprintf('fold.(ROIs{ri}).test.cc.%s = fold.(ROIs{ri}).test.cc.%s(cidx);',fs{fi},fs{fi}));
        end
        fs = fields(fold.(ROIs{ri}).test.time);
        for fi = 1:length(fs)
            eval(sprintf('fold.(ROIs{ri}).test.time.%s = fold.(ROIs{ri}).test.time.%s(tidx);',fs{fi},fs{fi}));
        end
        for fi = 1:length(fs)
            eval(sprintf('fold.(ROIs{ri}).train.time.%s = fold.(ROIs{ri}).train.time.%s(tidxtrain);',fs{fi},fs{fi}));
        end
    end
    
    dd.(sprintf('fold%i',cfold)) = fold; cfold=cfold+1;
end

dd.nfolds = cfold-1;
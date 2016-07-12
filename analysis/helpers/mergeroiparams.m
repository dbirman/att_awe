function params = mergeroiparams(roiparams,rois)

fs = fields(roiparams{1});
for fi = 1:length(fs)
    dat.(fs{fi}) = [];
end
for ri = 1:length(rois)
    roi = rois(ri);
    for fi = 1:length(fs)
        dat.(fs{fi}) = [dat.(fs{fi}) roiparams{roi}.(fs{fi})];
    end
end

for fi = 1:length(fs)
    params.(fs{fi}) = mean(dat.(fs{fi}));
end
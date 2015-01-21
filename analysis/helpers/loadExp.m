function expHolder = loadExp

expHolder = cell(length(files));
for fi = 1:length(files)
    fixFile = sprintf('%s/%s',analysis.datFolder,['f' files(fi).name]);
    load(fixFile);
    expHolder{fi} = exp;
end



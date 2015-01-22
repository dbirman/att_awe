function expHolder = loadExp(files)

global analysis

expHolder = cell(length(files));
for fi = 1:length(files)
    fixFile = sprintf('%s/%s',analysis.datFolder,['f' files(fi).name]);
    load(fixFile);
    expHolder{fi} = e;
end



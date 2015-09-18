function expHolder = loadExp(files)

global analysis

expHolder = cell(1,length(files));
for fi = 1:length(files)
    fixFile = sprintf('%s/%s',analysis.datFolder,['f' files(fi).name]);
    load(fixFile);
    expHolder{fi} = e;
end



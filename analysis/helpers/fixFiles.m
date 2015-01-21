function fixFiles(files)

global analysis

for fi = 1:length(files)
    curFile = sprintf('%s/%s',analysis.datFolder,files(fi).name);
    fixFile = sprintf('%s/%s',analysis.datFolder,['f' files(fi).name]);
    if isfile(fixFile)
        disp(sprintf('(noisecon_analysis) Found fixed file for %s, ignoring...',curFile));
    else
        success = genExp(curFile,fixFile);
        if success
            disp(sprintf('(noisecon_analysis) Output in: %s',fixFile));
        else
            disp(sprintf('(noisecon_analysis) File %s was not fixed!!!!',curFile));
        end
    end
end
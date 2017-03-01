function cc_copyFiles(sid)
% Input sid as 's300'

volDir = '/Volumes/data/cohcon';
if ~isdir(volDir)
    volDir = input('Input dubonnet directory: ','s');
end
behavTop = '~/data/cohcon';
disp('##############');
disp(sprintf('## SUBJ: %s',sid));
disp('##############');

behavDir = fullfile(behavTop,sid);
if ~isdir(behavDir)
    mkdir(behavDir);
end

mfiles = dir(sprintf('%s/%s/17*.mat',volDir,sid));
efiles = dir(sprintf('%s/%s/17*.edf',volDir,sid));

for mi = 1:length(mfiles)
    mfile = mfiles(mi);

    if ~isfile(fullfile(behavDir,mfile.name))
        disp(sprintf('Copying %s to %s',fullfile(volDir,sid,mfile.name),fullfile(behavDir,mfile.name)));
        copyfile(fullfile(volDir,sid,mfile.name),fullfile(behavDir,mfile.name));
    else
        disp(sprintf('Found file %s',fullfile(behavDir,mfile.name)));
    end
end
for ei = 1:length(efiles)
    efile = efiles(ei);

    if ~isfile(fullfile(behavDir,efile.name))
        disp(sprintf('Copying %s to %s',fullfile(volDir,sid,efile.name),fullfile(behavDir,efile.name)));
        copyfile(fullfile(volDir,sid,efile.name),fullfile(behavDir,efile.name));
    else
        disp(sprintf('Found file %s',fullfile(behavDir,efile.name)));
    end
end
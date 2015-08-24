function adFile = saveAllData( sid, allData )
%SAVEALLDATA Summary of this function goes here
%   Detailed explanation goes here

if ~isdir('~/proj/att_awe/analysis/data/')
    mkdir('~/proj/att_awe/analysis/data/');
end


adFile = fullfile('~/proj/att_awe/analysis/data/',sprintf('%s_allData.mat',sid));

if isfile(adFile)
    if ~isdir('~/proj/att_awe/analysis/data/adBackups')
        mkdir('~/proj/att_awe/analysis/data/adBackups');
    end
    copyfile(adFile,sprintf('~/proj/att_awe/analysis/data/adBackups/%s_%s_allData.mat',datestr(date,'YYMMDD'),sid));
end

disp('(allData) Saving...');
save(adFile,'allData');

end


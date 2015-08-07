function adFile = saveAllData( allData )
%SAVEALLDATA Summary of this function goes here
%   Detailed explanation goes here


adFile = '~/proj/att_awe/analysis/allData.mat';

if isfile(adFile)
    if ~isdir('~/proj/att_awe/analysis/adBackups')
        mkdir('~/proj/att_awe/analysis/adBackups');
    end
    copyfile(adFile,sprintf('~/proj/att_awe/analysis/adBackups/%s_allData.mat',datestr(date,'YYMMDD')));
end

disp('(allData) Saving...');
save(adFile,'allData');

end


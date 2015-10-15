function adFile = saveAllData( sid, allData )
%SAVEALLDATA Summary of this function goes here
%   Detailed explanation goes here

if ~isdir('~/data/cohcon/analysis/')
    mkdir('~/data/cohcon/analysis/');
end


adFile = fullfile('~/data/cohcon/analysis/',sprintf('%s_allData.mat',sid));


if isfile(adFile)
    if ~isdir('~/data/cohcon/analysis/adBackups')
        mkdir('~/data/cohcon/analysis/adBackups');
    end
    copyfile(adFile,sprintf('~/data/cohcon/analysis/adBackups/%s_%s_allData.mat',datestr(date,'YYMMDD'),sid));
end

disp('(allData) Saving...');
save(adFile,'allData');

end


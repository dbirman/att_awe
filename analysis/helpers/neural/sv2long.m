function [con, coh, task, sv] = sv2long ( sV, sN, tsV, tsN )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

con = []; coh = []; task=[]; sv= []; 

tasks = zeros(1,length(tsV));
for ti = 1:length(tsV)
    tasks(ti) = mrStr2num(tsN{ti}(end:end));
end
for i = 1:length(sV)
    sn = sN{i};
    csvs = sV{i};
    for ci = 1:length(csvs)
        csv = csvs(ci);
        sv = [sv csv];
        conPos = strfind(sn,'Con');
        cohPos = strfind(sn,'Coh');

        con = [con cc_getStringNum(sn,conPos)];
        coh = [coh cc_getStringNum(sn,cohPos)];

        % now find task
        for ti = 1:length(tsV)
            tsv = tsV{ti};
            if ~isempty(find(tsv==csv,1))
                task = [task tasks(ti)];
            end
        end
        if length(task) ~= length(con)
            disp('(cc_simplifySCM) Something is wrong with the stimVols for this run.');
            keyboard
        end
    end
end


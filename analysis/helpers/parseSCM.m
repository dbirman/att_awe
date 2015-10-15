
function [con, coh, task, newsV] = parseSCM(sN,sV,tsN, tsV,prefixs)

newsV = []; con = []; coh = []; task=[];

tasks = zeros(1,length(tsV));
for ti = 1:length(tsV)
    tasks(ti) = mrStr2num(tsN{ti}(end:end));
end
for i = 1:length(sV)
    if  ~isempty(sV{i})
        sn = sN{i};
        sv = sV{i};
        for si = 1:length(sv)
            csv = sv(si);
            newsV = [newsV csv];
            conPos = strfind(sn,prefixs{1});
            cohPos = strfind(sn,prefixs{2});
            
            con = [con cc_getStringNum(sn,conPos)];
            coh = [coh cc_getStringNum(sn,cohPos)];
            
            % now find task
            for ti = 1:length(tsV)
                tsv = tsV{ti};
                if ~isempty(find(tsv==csv,1))
                    task = [task tasks(ti)];
                end
            end
        end
    end
    if length(task) ~= length(con)
        disp('(cc_simplifySCM) Something is wrong with the stimVols for this run.');
        keyboard
    end
end

end
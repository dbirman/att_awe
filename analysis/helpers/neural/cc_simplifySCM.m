function nSCM = cc_simplifySCM( SCM, conValues, cohValues, nearest)

% We're going to modify the SCM for each folder in SCM.f## so that they
% only have values from conValues and cohValues.
% nearest = 0: Just those values specifically
% nearest = 1: Use the nereast neighbors

nSCM = struct;

fList = fields(SCM);
stims = {'lStim','rStim'};
prefixes = {{'rCon','rCoh'},{'lCon','lCoh'}};
pattern = '%s=%0.2f and %s=%0.2f and task=%i';
taskValues = [1 2];

if isempty(conValues)
    disp('(simplifySCM) Cannot simplify conValues: using all available values');
        
    % we need to build up a list of all THE VALUES
    for fi = 1:length(fList)
        folder = fList{fi};
        if folder(1)=='f'
            for si = 1:length(stims)
                stim = stims{si};
                prefixs = prefixes{si};

                sN = SCM.(folder).(stim).stimNames;
                sV = SCM.(folder).(stim).stimVol;
                tsN = SCM.(folder).(stim).taskNames;
                tsV = SCM.(folder).(stim).taskSV;

                % re-build into separate con/coh structures
                [con, ~, ~, ~] = parseSCM(sN,sV,tsN, tsV, prefixs);
                conValues = unique([conValues con]);
            end
        end
    end
end


if isempty(cohValues)
    disp('(simplifySCM) Cannot simplify cohValues: using all available values');
        
    % we need to build up a list of all THE VALUES
    for fi = 1:length(fList)
        folder = fList{fi};
        if folder(1)=='f'
            for si = 1:length(stims)
                stim = stims{si};
                prefixs = prefixes{si};

                sN = SCM.(folder).(stim).stimNames;
                sV = SCM.(folder).(stim).stimVol;
                tsN = SCM.(folder).(stim).taskNames;
                tsV = SCM.(folder).(stim).taskSV;

                % re-build into separate con/coh structures
                [~, coh, ~, ~] = parseSCM(sN,sV,tsN, tsV, prefixs);
                cohValues = unique([cohValues coh]);
            end
        end
    end
end

for fi = 1:length(fList)
    folder = fList{fi};
    disp(sprintf('Simplifying SCM for folder %s',folder));
    
    if folder(1)=='f'
        for si = 1:length(stims)
            stim = stims{si};
            prefixs = prefixes{si};

            sN = SCM.(folder).(stim).stimNames;
            sV = SCM.(folder).(stim).stimVol;
            tsN = SCM.(folder).(stim).taskNames;
            tsV = SCM.(folder).(stim).taskSV;

            % re-build into separate con/coh structures
            [con, coh, task, newsV] = parseSCM(sN,sV,tsN, tsV, prefixs);

            if ~isempty(conValues) && ~isempty(cohValues)
                % remove values we didn't ask for
                if nearest
                    % just move everything to the nearest values
                    con = nearestValue(con, conValues);
                    coh = nearestValue(coh, cohValues);
                else
                    % remove everything except if it's in conValues or cohValues
                    for ci = 1:length(con)
                        if ~any(conValues==con(ci)) || ~any(cohValues==coh(ci))
                            newsV(ci) = [];
                            con(ci) = [];
                            coh(ci) = [];
                        end
                    end
                end
            end

            % re-print stimNames and stimVol
            stimNames = {}; finalSV = {};
            for ccon = conValues
                for ccoh = cohValues
                    for ctask = taskValues
                        stimNames{end+1} = sprintf(pattern,prefixs{1},ccon,prefixs{2},ccoh,ctask);
                        % now we find corresponding stimVols
                        conP = con==ccon; cohP = coh==ccoh; taskP = task==ctask;
                        P = conP.*cohP.*taskP; % spots where we overlap
                        stimVolP = newsV(logical(P));
                        if isempty(stimVolP)
                            finalSV{end+1} = [];
                        else
                            finalSV{end+1} = sort(unique(stimVolP));
                        end
                    end
                end
            end

            % add to name
            nSCM.(folder).(stim).stimVol = finalSV;
            nSCM.(folder).(stim).stimNames = stimNames;
        end
    end
end

end

function [con, coh, task, newsV] = parseSCM(sN,sV,tsN, tsV,prefixs)

newsV = []; con = []; coh = []; task=[];

tasks = zeros(1,length(tsV));
for ti = 1:length(tsV)
    tasks(ti) = mrStr2num(tsN{ti}(end:end));
end
for i = 1:length(sV)
    if ~isempty(sV{i})
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

function val = nearestValue(val,values)

if length(val)>1
    for i = 1:length(val)
        val(i) = nearestValue(val(i),values);
    end
    return
end

if val > 1 || val < 0
    disp('only works on percentages');
    keyboard
end

if val==0, return; end

high = find(values >= val,1);
low = find(values <= val,1);

if isempty(low)
    val = values(1);
    return
end
if isempty(high)
    val = values(end);
    return
end

if abs(values(high)-val) < abs(values(low)-val)
    val = values(high);
else
    val = values(low);
end

end
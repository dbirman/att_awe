
function [side, rNum] = parseROI(roi,shortList)

rNum = -1;

sname = roi(1);


if strcmp(sname,'l')
    side = 1;
else
    side = 2;
end

found = [];
for i = 1:length(shortList)
    pos = strfind(roi,shortList{i});
    if ~isempty(pos)
        found(i) = pos;
    else
        found(i) = -1;
    end
end
rNum = find(found>0,1);
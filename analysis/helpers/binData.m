function [ Ybinned ] = binData( Y,X,bins )

Ybinned = {};
for i = 1:length(bins)-1
    low = bins(i);
    high = bins(i+1);
    glow = X>=low;
    shigh = X<high;
    Ybinned{end+1} = Y(find(logical(glow.*shigh)));

end


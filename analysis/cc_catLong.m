function fullData = cc_catLong(fdatas)

fullData = [];

for i = 1:length(fdatas)
    fullData = [fullData ; fdatas{i}];
end
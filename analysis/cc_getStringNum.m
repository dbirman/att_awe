

function num = cc_getStringNum(str,j)
k = j+1;
while j <= length(str) && k <= length(str) && isempty(str2num(str(j:k)))
    j = j+1; k = k+1;
end
while k <= length(str) && ~isempty(str2num(str(j:k)))
    k = k+1;
end

num = str2num(str(j:k-1));
end
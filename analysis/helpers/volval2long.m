function data = volval2long( vol, val )
%VOLVAL2LONG Convert volume and value data to long form
%
% data = volval2long(vol,val);

data = zeros(10000,2); count = 1;
for ci = 1:length(val)
    for vi = 1:length(vol{ci})
        data(count,:) = [vol{ci}(vi) val(ci)];
        count = count + 1;
    end
end
data = data(1:count-1,:);


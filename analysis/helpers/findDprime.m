function deltaX = findDprime( thresh, x, resp, poisson )
% Tries to find the deltaX value such that thresh+deltaX = resp(thresh)+1

deltaX = zeros(size(thresh));
if length(thresh)>1
    for i = 1:length(thresh)
        deltaX(i) = findDprime(thresh(i),x,resp,poisson);
    end
    return
end

deltaX = [];

if poisson
    inc = sqrt(abs(resp(find(x>=thresh,1))));
else
    inc = 1/(sqrt(2)*norminv(0.815));
end

while isempty(deltaX);
    % thresh is 1 value
    deltaX = x(find(resp>=(inc+resp(find(x>=thresh,1))),1))-thresh;
    if isempty(deltaX)
        deltaX = NaN;
        break
%         inc = inc*0.95;
%         if inc<0.001
%             deltaX = NaN;
%             break
%         end
    end
end

% deltaX = deltaX/inc;
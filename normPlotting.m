function [plotting plottingsd] = normPlotting(plotting)

plotcell = plotting;
plotting = zeros(2,3,2,2);
plottingsd = zeros(2,3,2,2);

for i = 1:3
    for j = 1:2
        dat = plotcell{1,i,1,j};
        dat = dat(logical([dat<.55].*[dat>0]));
        norms(i,j) = mean(dat);
    end
end

% reorganize plotting
for i = 1:2
    for j = 1:3
        for k = 1:2
            for l = 1:2
                dat = plotcell{i,j,k,l};
                dat = dat(logical([dat<.55].*[dat>0]));
                if norms
                    dat = dat / norms(j,l);
                end
                plotting(i,j,k,l) = mean(dat);
                plottingsd(i,j,k,l) = std(dat)/sqrt(length(dat));
            end
        end
    end
end
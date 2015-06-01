%% Script to test naka-rushton output...

sides = {'left'};
ROIs = {'l_v1','l_hmt'};

cueds = {'uncued','cued'};
values = {'contrast','coherence'};

%%
linetypes = {'-','--'};

for si = 1:length(sides)
    for ri = 1:length(ROIs)
        cData = allData.(sides{si}).(ROIs{ri});
        
        legVals = {};
        figure, hold on
        for vi = 1:length(values)
            for cui = 1:length(cueds)
                legVals{end+1} = sprintf('%s.%s',values{vi},cueds{cui});
                
                c = cData.(values{vi}).(cueds{cui}).c;
                r = cData.(values{vi}).(cueds{cui}).r;
                n = cData.(values{vi}).(cueds{cui}).n;
                
                c = c(r>0);
                n = n(r>0);
                r = r(r>0);
                
                curCol = rand(1,3);
                
                plot(c,r,'*','color',curCol);
                
                try
%                     fit = fitNakaRushton(c,r);
                fit = fitNakaRushtonWeighted(c,r,n);
                k = fit.Rmax;
                n = fit.n;
                cn = c;
                c50 = fit.c50;
                
                y = k*(c.^n./(cn.^n+c50.^n))+fit.offset;
                
                plot(c,y,linetypes{vi},'color',curCol);
                legVals{end+1} = sprintf('%s.%s.fit',values{vi},cueds{cui});
                catch
                end
            end
        end
        legend(legVals);
        title(sprintf('%s.%s',sides{si},ROIs{ri}));
    end
end
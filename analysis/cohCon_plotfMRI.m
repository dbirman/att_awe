%% Script to test naka-rushton output...

ROIs = {'v1','hmt'};

cueds = {'uncued','cued'};
values = {'contrast','coherence'};

%% Plot amp responses

for si = 1:length(sides)
    for ri = 1:length(ROIs)
        cData = allData.(sides{si}).(ROIs{ri});
        
        legVals = {};
        for vi = 1:length(values)
            for cui = 1:length(cueds)   
                figure, hold on             
                legVals{end+1} = sprintf('%s.%s',values{vi},cueds{cui});
                
                amps = cData.(values{vi}).(cueds{cui}).a;
                for ai = 1:length(amps)
                    plot(cData.(values{vi}).(cueds{cui}).time,cData.(values{vi}).(cueds{cui}).canon.*amps(ai));
                end
            end
        end
        
    end
end

%%
linetypes = {'-','--'};

figure
for ri = 1:length(ROIs)
    cData = allData.(ROIs{ri});

    for vi = 1:length(values)
        legVals = {};
        subplot(length(ROIs),length(values),(vi-1)*length(ROIs)+ri), hold on
        for cui = 1:length(cueds)
            legVals{end+1} = sprintf('%s.%s',values{vi},cueds{cui});

            c = cData.(values{vi}).(cueds{cui}).i;
            r = cData.(values{vi}).(cueds{cui}).a;
            n = cData.(values{vi}).(cueds{cui}).N;
            se = cData.(values{vi}).(cueds{cui}).ase;
            c = c(n>40);
            r = r(n>40);
            se = se(n>40);
            n = n(n>40);
            
            ci = real(se) ./sqrt(n).*1.96;


            curCol = rand(1,3);

            errorbar(c,r,se,'-*','color',curCol);

        end
        legend(legVals);
        title(sprintf('%s',ROIs{ri}));
    end
end
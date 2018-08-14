%% Plot beta values
betas = zeros(length(riopts),3);
for ri = 1:8
    ropt = riopts(ri);
    
    y = cc_(ri,:);%y = fliplr(y);
    con = avg_decon.conidx;
    coh = avg_decon.cohidx;
    con = con-min(con);
    
    betas(ri,:) = [ones(size(con')) con' coh']\y'; % additive model
 
    conidx = avg_decon.conidx;
    cohidx = avg_decon.cohidx;
    ucon = unique(conidx);
    ucoh = unique(cohidx);
    
    h = figure;
    
    temp = betas(ri,:);
    
    cmap = brewermap(11,'PuOr');
    mapcon = [5 4 3 2];
    mapcoh = [7 8 9 10 11];
    % two subplots, one for contrast responses (across coherence) one for
    % coherence responses (across contrast)
    subplot(211); hold on
    title('Effect of increasing contrast at different coherences');
    % contrast across coherence
    legs = {};
    for ui = 1:length(ucoh)
        % get the values corresponding to the lowest level of coherence
        vals = cohidx==ucoh(ui);
        plot(ucon*100,y(vals)','Color',cmap(mapcoh(ui),:));
%         errbar(ucon*100,y(vals)',ccs(vals)','Color',cmap(mapcoh(ui),:));
        legs{end+1} = sprintf('Coherence +%i%%',ucoh(ui)*100);
    end
%     legend(legs)
    xlabel('Contrast (%)');
    ylabel('Response (% signal change)');
    
    subplot(212); hold on
    title('Effect of increasing coherence at different contrasts');
    % contrast across coherence
    legs = {};
    for ui = 1:length(ucon)
        % get the values corresponding to the lowest level of coherence
        vals = conidx==ucon(ui);
        plot(ucoh*100,y(vals)','Color',cmap(mapcon(ui),:));
%         errbar(ucoh*100,y(vals)',ccs(vals)','Color',cmap(mapcon(ui),:));
        legs{end+1} = sprintf('Contrast +%i%%',ucon(ui)*100);
    end
%     legend(legs)
    xlabel('Coherence (%)');
    ylabel('Response (% signal change)');
    
%     savepdf(h,fullfile(datafolder,sprintf('avg_effect_%s.pdf',ROIs{ropt})));
end

coherences = [0.15 0.3 0.45 0.6];
thresh_vals = .05:.05:.65;
contrasts = [0.325 0.4 0.55 0.85];
% errbar(coherences,control(1,:),controls(1,:),'Color',map(6,:));
% errbar(contrasts,control(2,:),controls(2,:),'Color',map(2,:));

h1 = plot(coherences,squeeze(control(ai,1,:)),'o','MarkerSize',10,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(6,:));
h5 = plot(contrasts,squeeze(control(ai,2,:)),'o','MarkerSize',10,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',map(2,:));

try
if ii==1
    mylegend({'Coherence','Contrast'},{{'o','Color',map(6,:),'MarkerFaceColor',map(6,:)},{'o','Color',map(2,:),'MarkerFaceColor',map(2,:)}},'boxoff=1');
    l = legend(gca,'boxoff');
    set(l,'Color','none');
end

if ii==4
    xlabel('Base strength (%)');
    ylabel('JND (%)');
end
catch
end
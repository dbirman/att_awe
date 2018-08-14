

%% Plotting

for pi = 1:2
    for li = 1:3
        errbar(contrasts(pi)+li*.005,thresh_(2,pi,li),squeeze(thresh(2,2,pi,li))-thresh_(2,pi,li),'-','Color',cmap(1,:));
        plot(contrasts(pi)+li*.005,thresh_(2,pi,li),'o','MarkerFaceColor',cmap(1,:),'MarkerSize',li*2,'MarkerEdgeColor','white');
        
        errbar(coherences(pi)+li*.005,thresh_(1,pi,li),squeeze(thresh(2,1,pi,li))-thresh_(1,pi,li),'-','Color',cmap(3,:));
        plot(coherences(pi)+li*.005,thresh_(1,pi,li),'o','MarkerFaceColor',cmap(3,:),'MarkerSize',li*2,'MarkerEdgeColor','white');
    end
end
p(1) = plot(-1,-1,'ok','MarkerFaceColor','black','MarkerSize',2);
p(2) = plot(-1,-1,'ok','MarkerFaceColor','black','MarkerSize',4);
p(3) = plot(-1,-1,'ok','MarkerFaceColor','black','MarkerSize',6);
legend(p,{'250 ms','500 ms','1000 ms'});
axis([0 0.75 0 0.4]);
set(gca,'YTick',[0 .1 .2 .3],'YTickLabel',[0 10 20 30]);
set(gca,'XTick',[0 0.25 0.5 0.75],'XTicklabel',{'0%','25%','50%','75%'});
xlabel('\Delta stimulus (%)');
ylabel('JND (%)');

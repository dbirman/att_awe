function [f, data] = cc_stateTracePlot( sid )

allData = loadAllData(sid);

data = struct;
x = 0:.0001:1;

atts = {'att_con','att_coh'};
fits = {'cohfit','confit'};

data.att_con.cohfit = allData.behav.caf{1}.fit;
data.att_con.confit = allData.behav.maf{2}.fit;
data.att_coh.cohfit = allData.behav.maf{1}.fit;
data.att_coh.confit = allData.behav.caf{2}.fit;

%%
for ai = 1:length(atts)
    att = atts{ai};
    for fi = 1:length(fits)
        fit = fits{fi};
        data.(att).(fit).st_x = x;
        data.(att).(fit).st_y = weibull(x,data.(att).(fit).fitparams);
    end
end

%%

for ai = 1:length(atts)
    att = atts{ai};
    data.(att).coh = [];
    data.(att).con = [];
    data.(att).p = [];
    for p = .5:.001:1
        valcoh = data.(att).cohfit.st_x(find(data.(att).cohfit.st_y>p,1));
        valcon = data.(att).confit.st_x(find(data.(att).confit.st_y>p,1));
        if ~isempty(valcoh) && ~isempty(valcon)
            data.(att).p(end+1) = p;
            data.(att).coh(end+1) = valcoh;
            data.(att).con(end+1) = valcon;
        end
    end
end

%% corr

for ai = 1:length(atts)
    att = atts{ai};
    
    data.(att).r = corr(data.(att).coh',data.(att).con','type','Pearson');
    data.(att).sr = corr(data.(att).coh',data.(att).con','type','Spearman');
end
    
%%

f = figure();
hold on
plot(data.att_con.con,data.att_con.coh,'*r');
plot(data.att_coh.con,data.att_coh.coh,'*b');
legend({'Attending Contrast','Attending Coherence'});

% add the lines to connect equal p dots
for p = .5:.025:1
    i = find(data.att_con.p==p);
    j = find(data.att_coh.p==p);
    if ~isempty(j) && ~isempty(i)
        plot([data.att_con.con(i) data.att_coh.con(j)],[data.att_con.coh(i) data.att_coh.coh(j)],'y');
        text(abs(data.att_coh.con(j)-data.att_con.con(i))/2+min(data.att_con.con(i),data.att_coh.con(j)),abs(data.att_con.coh(i)-data.att_coh.coh(j))/2+min(data.att_con.coh(i),data.att_coh.coh(j)),num2str(p));
    elseif ~isempty(i)
        % we can add for att_con
        text(data.att_con.con(i),data.att_con.coh(i)-.01,num2str(p));
    elseif ~isempty(j)
        % we can add for att_coh
        text(data.att_coh.con(j),data.att_coh.coh(j)-.01,num2str(p));
    end
end

x = 0:.001:1;
maf2t.y = weibull(x,allData.s304.behav.maf{2}.fit.fitparams);
maf2t.t = x(find(maf2t.y>=.816,1));
caf1t.y = weibull(x,allData.s304.behav.caf{1}.fit.fitparams);
caf1t.t = x(find(caf1t.y>=.816,1));
caf2t.y = weibull(x,allData.s304.behav.caf{2}.fit.fitparams);
caf2t.t = x(find(caf2t.y>=.816,1));
maf1t.y = weibull(x,allData.s304.behav.maf{1}.fit.fitparams);
maf1t.t = x(find(maf1t.y>=.816,1));

plot(maf2t.t,caf1t.t,'*k');
plot(caf2t.t,maf1t.t,'*k');

% compute spearman's r

allcoh = [data.att_con.coh data.att_coh.coh];
allcon = [data.att_con.con data.att_coh.con];

[allcon, i] = sort(allcon);
allcoh = allcoh(i);

% plot(allcon,allcoh)

r = corr(allcon',allcoh','type','Spearman');
r

title('State-Trace Plot: Points Correspond to a %Correct Rate');
xlabel('Contrast (Delta Intensity)');
ylabel('Coherence (Delta Intensity)');

set(gca,'FontSize',14);

%% Testing
% 
% x = 0:.0001:1;
% wfit = allData.s304.behav.caf{1}.fit.fitparams;
% % wfit(3) = 0; wfit(4) = .5;
% y = weibull(x,wfit);
% 
% 
% figure
% 
% plot(x,y);
% 
% p = y(find(x>=wfit(1),1))

%%
t = 0:pi/60:2*pi;

xs = randi(500,1,500);
ys = randi(500,1,500);

con = rand(1,500)<0.5;

h= figure;
iptsetpref('ImshowBorder','tight');
for ti = 1:length(t)
    ct = t(ti);
    ccoh = 0.60;%(1 + cos(ct/2))/2;
    ccon = 1;%(1 + cos(ct))/2;
    cdir = 2;%cos(ct/5);
    
    clf; hold on
    blk = [0.5 0.5 0.5]-0.5*ccon;
    wht = [0.5 0.5 0.5]+0.5*ccon;
    plot(xs(con),ys(con),'ok','MarkerFaceColor',blk,'MarkerEdgeColor',blk);
    plot(xs(~con),ys(~con),'ow','MarkerFaceColor',wht,'MarkerEdgeColor',wht);
    set(gca,'Color',[0.5 0.5 0.5]);
    pause(.001);
    
%     axis square
    
    % update location
    idx = rand(1,500)<ccoh;
    xs(idx) = xs(idx)+cdir;
    rt = rand(1,sum(~idx))*2*pi;
    xs(~idx) = xs(~idx) + cdir*cos(rt);
    ys(~idx) = ys(~idx) + cdir*sin(rt);
    
    ys(ys<0)=ys(ys<0)+500;
    ys(ys>500)=ys(ys>500)-500;
    xs(xs<0)=xs(xs<0)+500;
    xs(xs>500)=xs(xs>500)-500;
    
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    
    pause(.01);
    F = getframe;
    [X,map] = frame2im(F);
    imwrite(X,fullfile(sprintf('~/proj/att_awe/talks/data_figures/coh60/cohcon_%03.0f.png',ti)));
end
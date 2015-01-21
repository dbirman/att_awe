function gen_perf( stimulus , plotting, peripheral)

for dual = 1:2
    if dual == 1
        stairtype = 'staircase';
        typeD = 'single';
    else
        stairtype = 'dualstaircase';
        typeD = 'dual';
    end
    for type = 1:2
        if type == 1
            dispText = 'Noise';
        else
            dispText = 'Contrast';
        end
        main.(dispText).(typeD) = mean(plotting(:,:,dual,type),2);
    end
end

% Main task performance
mainConFocalPerf = main.Contrast.single(1);
mainConDistPerf = main.Contrast.single(2);
mainNoiseFocalPerf = main.Noise.single(1);
mainNoiseDistPerf = main.Noise.single(2);

mainConFocalPerf = mainConFocalPerf / mainConFocalPerf;
mainConDistPerf = mainConDistPerf / mainConDistPerf;
mainNoiseFocalPerf = mainNoiseFocalPerf / mainNoiseFocalPerf;
mainNoiseDistPerf = mainNoiseDistPerf / mainNoiseDistPerf;

mainConFocalDualPerf = main.Contrast.dual(1);
mainConDistDualPerf = main.Contrast.dual(2);
mainNoiseFocalDualPerf = main.Noise.dual(1);
mainNoiseDistDualPerf = main.Noise.dual(2);

mainConFocalDualPerf = mainConFocalDualPerf / mainConFocalPerf;
mainConDistDualPerf = mainConDistDualPerf / mainConDistPerf;
mainNoiseFocalDualPerf = mainNoiseFocalDualPerf / mainNoiseFocalPerf;
mainNoiseDistDualPerf = mainNoiseDistDualPerf / mainNoiseDistPerf;

for i = 1:3
    cur = peripheral{i};
    cur = cur(logical([cur<.25].*[cur>0]));
    peripheral{i} = cur;
end

gPerf = mean(peripheral{1});
gNPerf = mean(peripheral{2});
gCPerf = mean(peripheral{3});

gPerf = gPerf / gPerf;
gNPerf = gNPerf / gPerf;
gCPerf = gCPerf / gCPerf;

% Plot
figure
hold on
title('Dual Task Performance');
% singles
plot(0,gPerf,'*r');
ylabel('Gender (SOA ms)');
xlabel('Contrast/Noise Performance (delta)');
plot(mainNoiseFocalPerf,0,'*g');
text(mainNoiseFocalPerf,.01,sprintf('%0.2f',mainNoiseFocalPerf/(1-mainNoiseFocalPerf)));
plot(mainNoiseDistPerf,0,'*c');
text(mainNoiseDistPerf,.01,sprintf('%0.2f',mainNoiseDistPerf/(1-mainNoiseDistPerf)));
plot(mainConFocalPerf,0,'*r');
plot(mainConDistPerf,0,'*m');
% duals
plot(mainNoiseFocalDualPerf,gNPerf,'*g');
text(mainNoiseFocalDualPerf,gNPerf+.01,sprintf('%0.2f',mainNoiseFocalDualPerf/(1-mainNoiseFocalDualPerf)));
plot(mainNoiseDistDualPerf,gNPerf,'*c');
text(mainNoiseDistDualPerf,gNPerf+.01,sprintf('%0.2f',mainNoiseDistDualPerf/(1-mainNoiseDistDualPerf)));
plot(mainConFocalDualPerf,gCPerf,'*r');
plot(mainConDistDualPerf,gCPerf,'*m');
% lines
plot(0:mainNoiseDistPerf/10:mainNoiseDistPerf,repmat(gPerf,1,11),'--r');
plot(0:.25:.5,repmat(.25,1,3),'--k');
plot(repmat(mainNoiseFocalPerf,1,11),0:gPerf/10:gPerf,'--g');
plot(repmat(mainNoiseDistPerf,1,11),0:gPerf/10:gPerf,'--c');
plot(repmat(mainConFocalPerf,1,11),0:gPerf/10:gPerf,'--r');
plot(repmat(mainConDistPerf,1,11),0:gPerf/10:gPerf,'--m');
plot(repmat(.5,1,3),0:.125:.25,'--k');
% legend
% legend('Gender','Focal Noise','Dist Noise','Focal Contrast','Dist Contrast');
% distance lines
xp = [mainConFocalPerf mainConFocalDualPerf];
yp = [gPerf gCPerf];
distCF = sqrt((xp(2)-xp(1))^2 + (yp(2)-yp(1))^2);
plot(xp,yp,'-k');
text(mean(xp),mean(yp),sprintf('Dist: %1.2f',distCF));

xp = [mainNoiseFocalPerf mainNoiseFocalDualPerf];
yp = [gPerf gNPerf];
distNF = sqrt((xp(2)-xp(1))^2 + (yp(2)-yp(1))^2);
plot(xp,yp,'-k');
text(mean(xp),mean(yp),sprintf('Dist: %1.2f',distNF));

xp = [mainConDistPerf mainConDistDualPerf];
yp = [gPerf gCPerf];
distCD = sqrt((xp(2)-xp(1))^2 + (yp(2)-yp(1))^2);
plot(xp,yp,'-k');
text(mean(xp),mean(yp),sprintf('Dist: %1.2f',distCD));

xp = [mainNoiseDistPerf mainNoiseDistDualPerf];
yp = [gPerf gNPerf];
distND = sqrt((xp(2)-xp(1))^2 + (yp(2)-yp(1))^2);
plot(xp,yp,'-k');
text(mean(xp),mean(yp),sprintf('Dist: %1.2f',distND));

disp(sprintf('Euclidian Distance for Contrast: %1.2f, for Noise: %1.2f',sqrt(distCF^2+distCD^2),sqrt(distNF^2+distND^2)));
disp('This distance measurement is innacurate :), noise should be on a log scale and euclidian distance isn''t really correct here');
 hold off


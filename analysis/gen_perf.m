function gen_perf( stimulus , plotting)

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

mainConFocalDualPerf = main.Contrast.dual(1);
mainConDistDualPerf = main.Contrast.dual(2);
mainNoiseFocalDualPerf = main.Noise.dual(1);
mainNoiseDistDualPerf = main.Noise.dual(2);

% Get the peripheral task performance
%%%% check for >2 task sets
if length(stimulus.p.dualstaircase{1}) > 1
    genderNoisePerf = doStaircase('threshold',stimulus.p.dualstaircase{1}(3:end),'type','weibull');
%     genderNoisePerf = doStaircase('threshold',stimulus.p.dualstaircase{1}(3:end));
else
    genderNoisePerf = doStaircase('threshold',stimulus.p.dualstaircase{1},'type','weibull');
%     genderNoisePerf = doStaircase('threshold',stimulus.p.dualstaircase{1});
end
if length(stimulus.p.dualstaircase{2}) > 1
    genderConPerf = doStaircase('threshold',stimulus.p.dualstaircase{2}(3:end),'type','weibull');
%     genderConPerf = doStaircase('threshold',stimulus.p.dualstaircase{2}(3:end));
else
    genderConPerf = doStaircase('threshold',stimulus.p.dualstaircase{2},'type','weibull');
%     genderConPerf = doStaircase('threshold',stimulus.p.dualstaircase{2});
end
% genderPerf = doStaircase('threshold',stimulus.p.staircase(2:end),'type','weibull');
genderPerf = doStaircase('threshold',stimulus.p.staircase(1:end-1),'type','weibull');
% gNPerf = genderNoisePerf.threshold;
% gCPerf = genderConPerf.threshold;
% gPerf = genderPerf.threshold;
gNPerf = genderNoisePerf.threshold;
gCPerf = genderConPerf.threshold;
gPerf = genderPerf.threshold;
% % % % % % % % % Normalize
% % % % % % % % gNPerf_N = gNPerf;
% % % % % % % % gCPerf_N = gCPerf;
% % % % % % % % gPerf = gPerf;
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


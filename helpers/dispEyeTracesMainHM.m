% % % % % % Analysis ideas:
% % % % % % 
% % % % % %  - Add the actual images at pedestals + threshold differences to the discrimination function image
% % % % % %  - Do normalisation, take 1/Normed, and then plot those to see the tradeoff
% % % % % %  - Get the estimated ?threshold? across every trial for both tasks during dual and plot them to see the ?tradeoff? curves. Plot both tradeoff curves next to each other
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % what are you currently working on:
% % % % % % 
% % % % % % getting eye tracking data to show properly? right now it displays sort of nicely, but it needs to be concatenated across runs. Realistically it needs to be done in R? Matlab just sucks!
%%%%%%%%%%%%%%%%%%%%%%%%%%
%    displayEyeTraces    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispEyeTracesMainHM(e)

% Heatmap version of dispEyeTraces

%%
main = e{1};
eye = e{4};

% Stimulus were presented in radial angles 3.8 deg away (at 45 degrees) and
% with diameter 5.5. We will calculate a square
stimSize = 5.5;
stimPos = 3.8;
H = sqrt(stimPos^2+stimPos^2);
cx = H*cos(deg2rad(45));
cy = H*sin(deg2rad(45));

segNames = {'ITI','Cue','Disp1','Resp1','Disp2','Resp2','RespMain'};
figure
set(gcf, 'units', 'centimeters', 'pos', [25 25 28 40])
meanX =[];
meanY=[];
for loc = 1:4
    for s = 1:7
        subplot(7,4,loc+(s-1)*4)
        hold on
        for t = 1:length(main.randVars.targetLoc)
            if loc==main.randVars.targetLoc(t)
                cT = main.trials(t);
                xData = eye.xPos(t,:);
                meanX = [meanX nanmean(xData)];
                yData = eye.yPos(t,:);
                meanY = [meanY nanmean(yData)];
                times = cT.segtime;
                times = times - times(1);
                % Looking at a specific segment
                lowBound = times(s);
                try
                    upBound = times(s+1);
                catch
                    upBound = inf;
                end
                pos = find([eye.time > lowBound] .* [eye.time < upBound]);
                plot(xData(pos),yData(pos),'*');
            end
        end
        title(segNames{s});
        axis([mean(meanX)-9 mean(meanX)+9 mean(meanY)-9 mean(meanY)+9])
        square(cx+mean(meanX),cy+mean(meanY),stimSize);
        square(mean(meanX)-cx,mean(meanY)-cy,stimSize);
        square(mean(meanX)-cx,cy+mean(meanY),stimSize);
        square(cx+mean(meanX),mean(meanY)-cy,stimSize);
        
        hold off
    end
end

function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp);

function square(x,y,e)
% x and y are the center
% e is the edge length

x0 = x - e/2; % top left
y0 = y - e/2; % top left
x1 = x + e/2;
y1 = y + e/2;
xs = 0:.5:e;
hold on
plot(xs+x0,repmat(y0,1,length(xs)),'r'); % top line
plot(xs+x0,repmat(y1,1,length(xs)),'r'); % bottom line
plot(repmat(x0,1,length(xs)),y0+xs,'r');
plot(repmat(x1,1,length(xs)),y0+xs,'r');
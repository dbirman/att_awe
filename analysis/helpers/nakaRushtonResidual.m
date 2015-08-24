
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    nakaRushtonResidual    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual = nakaRushtonResidual(params,c,r,parameterIndex,m)

dispFit = 1;
 
% calculate residual for each crf
for i = 1:size(c,1)
  thisC = c(i,:);
  if ~m.fixedN
    % decode parameters
    Rmax = params(parameterIndex(i,1));
    c50 = params(parameterIndex(i,2));
    n = params(parameterIndex(i,3));
    offset = params(parameterIndex(i,4));
    % calculate naka-rushton, note that we use the log of c
    fitR(i,:) = nakaRushton(thisC,[Rmax c50 n offset]); 
  else %fixedN
    % decode parameters
    Rmax = params(parameterIndex(i,1));
    c50 = params(parameterIndex(i,2));
    n = m.fixedN;
    offset = params(parameterIndex(i,3));
    % calculate naka-rushton, note that we use the log of c
    fitR(i,:) = nakaRushton(thisC,[Rmax c50 n offset]); 
  end
end

% display fit if called for
if dispFit
  f = smartfig('selectionModel_nakaRushtonResidual','reuse');  
  clf;
  for i = 1:size(c,1)
    subplot(1,size(c,1),i)
    semilogx(c(i,:),r(i,:),'ko');
    hold on
    semilogx(c(i,:),fitR(i,:),'k-')
    if ~m.fixedN
      titleStr = sprintf('Rmax: %0.3f c50: %0.2f n: %0.3f\n',params(parameterIndex(i,1)),params(parameterIndex(i,2)),params(parameterIndex(i,3)));
      title(sprintf('%s offset: %f',titleStr,params(parameterIndex(i,4))));
    else
      titleStr = sprintf('Rmax: %0.3f c50: %0.2f n: %0.3f\n',params(parameterIndex(i,1)),params(parameterIndex(i,2)),m.fixedN);
      title(sprintf('%s offset: %f',titleStr,params(parameterIndex(i,3))));
    end
  end
  makeEqualYaxis(1,size(c,1));
  drawnow
end

residual = r-fitR;
residual = residual(:);

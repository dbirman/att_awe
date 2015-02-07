function stimulus = initFadeGratings(stimulus,myscreen)

if ~isfield(stimulus,'gratingsInitialized')

    % calculate all the orientations we need
    stimulus.grating.orientations = 0:360;

    % calculate all the phases we are going to compute
    stimulus.grating.centerPhase = 0;
    stimulus.grating.nPhases = 36;
    stimulus.grating.phases = [0:2*pi/(stimulus.grating.nPhases):2*pi];
    stimulus.grating.phases = stimulus.grating.phases(1:end-1);
    nPhases = length(stimulus.grating.phases);

    % this gives the phase order to go back and forth
    stimulus.grating.phaseIndex = [1:nPhases 1 nPhases:-1:2];
    stimulus.grating.phaseIndexLen = length(stimulus.grating.phaseIndex);

    % make all the 1D gratings. We compute all phases and all possible contrast values given the
    % range of indexes available to us. The 1st texture is gray the nth texture is full
    % contrast for the current gamma setting
    stimulus.colors.nDisplayContrasts = floor(stimulus.colors.rmed);
    for iPhase = 1:nPhases
      for iContrast = 0:stimulus.colors.nDisplayContrasts
        pDone = calcPercentDone(iPhase,nPhases,iContrast,stimulus.colors.nDisplayContrasts);
        disppercent(pDone);
        if myscreen.userHitEsc,mglClose;keyboard,end
        % get the phase
        thisPhase = (stimulus.grating.centerPhase+stimulus.grating.phases(iPhase))*180/pi;
        % make the grating
        thisGrating = round(iContrast*mglMakeGrating(stimulus.grating.width,nan,stimulus.grating.sf,0,thisPhase)+stimulus.colors.rmed);
        if min(thisGrating(:)) < 0
          keyboard
        end
        % create the texture
        stimulus.gratings(iContrast+1,iPhase,:,:) = thisGrating;
      end
    end
    disppercent(inf);
    stimulus.gratingsInitialized = 1;
end

for i = 1:size(stimulus.gratings,1)
    for j = 1:size(stimulus.gratings,2)
        stimulus.tex(i,j) = mglCreateTexture(squeeze(stimulus.gratings(i,j,:,:)));
    end
end
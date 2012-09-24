function p = GetParams(action,ROOT)
%  Constructs a parameters set for Botrylloides
  
  %% Setup up simulation parameters
 
  simParams.numTailSegs =      200; 
  simParams.numTrunkSegs =     50;
  simParams.numTailBeats =     2; %500;
  simParams.sampleRate =       100; %(1/tailBeat)
  
  
  %% Get kinematic and morphological parameters
  
  switch action
      % Parameters for ascidian larvae
      case 'botry'
          [finParams,meatParams,ocellParams,kinParams]...
            = GetBotry(simParams,ROOT);
        
       % Parameters for zebrafish larvae
       case 'danio'
          [finParams,meatParams,ocellParams,kinParams]...
            = GetDanio(simParams,ROOT);
  end
  
  
  %% Setup time information

  p.timeStart = 0;
  p.timeStop = kinParams.stg1.dur + max(finParams.s)/kinParams.stg1.waveSpeed ...
               + kinParams.stg2.beatPeriod + ...
               simParams.numTailBeats*kinParams.und.beatPeriod;
  
  
  %% Setup light

  p.lightT = linspace(p.timeStart,p.timeStop,...
      simParams.numTailBeats*simParams.sampleRate);
  p.light = [zeros(1,length(p.lightT)); 
     zeros(1,length(p.lightT)); -ones(1,length(p.lightT))];
         
  p.dlight = [zeros(3,1) diff(p.light,1,2)];
  
  
  %% Setup gravity
  
  p.gravity = [0; 0; -9.81];
  
  
  %% Setup fluid constants
  
  p.fluidDensity = 1.02e-2;    % qg/qm^3
  p.fluidkVisc = 1.047e0;       % qm^2/cs
  

  %% Calculate tail information
    
  % get kinematics
  [p.larvaTailT,p.larvaTailS,p.larvaTailRX,p.larvaTailRY] = ...
      GetKinematics(kinParams,simParams,max(finParams.s),p.timeStop);
  
  % make periodic version
  s = linspace(min(p.larvaTailS), max(p.larvaTailS), length(p.larvaTailS));
  t = linspace(min(p.larvaTailT), 3*max(p.larvaTailT), 3*length(p.larvaTailT));
  rx = [p.larvaTailRX; p.larvaTailRX; p.larvaTailRX];
  ry = [p.larvaTailRY; p.larvaTailRY; p.larvaTailRY];
  
  % spline fitting
  spRx = csaps({t, s}, rx, 0.999);
  spVx = fnder(spRx,[1 0]);
  spAx = fnder(spRx,[2 0]);
  
  spRy = csaps({t, s}, ry,.9999);
  spVy = fnder(spRy,[1 0]);
  spAy = fnder(spRy,[2 0]);
  
  % get domain of interrogation
  s = p.larvaTailS+10*eps;
  t = p.larvaTailT+max(p.larvaTailT)+10*eps;
 
  % get velocity and acceleration
  p.larvaTailUX = fnval(spVx,{t, s});
  p.larvaTailUY = fnval(spVy,{t, s});

  p.larvaTailAX = fnval(spAx,{t, s});
  p.larvaTailAY = fnval(spAy,{t, s});

  clear s t rx ry spRx spRy spVx spVy spAx spAy;
%  

  
  %% Copy parameters
  
  p.ocellParams = ocellParams;
  p.finParams = finParams;
  p.meatParams = meatParams;
  p.kinParams = kinParams;




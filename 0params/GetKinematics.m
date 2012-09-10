function [t,s,x,y] = getKinematics(kinParams,simParams,tailLength)
% This function returns the tail kinematics for every point
% along the tail for every instant in a tail beat
%  
% Input data structure
%  kinParams -> kLeft - scalar
%            -> kRight - scalar
%            -> mLeft - scalar
%            -> mRight - scalar
%            -> thetaAmp - scalar
%            -> waveSpeed - scalar
%            -> beatFreq - scalar
%
%  simParams -> numTailSegs - scalar
%            -> numTrunkSegs - scalar
%            -> numTailBeats - scalar
%            -> sampleRate - scalar
%
%  tailLength - scalar
%  
% Output data structure
%  t,s,x,y are matrices with the arc length 
%  in rows and the time in columns
%

  % Set oversample factor
  upsamp = 10;

  % Generate arc length, time grid
  [s,t] = meshgrid(linspace(0,tailLength,simParams.numTailSegs+1),...
           linspace(0,1/kinParams.beatFreq,upsamp*simParams.sampleRate));
  
  % Calculate ds
  ds = diff(s,1,2);
  ds = [ds(:,1) ds];
  
  % Calculate wave position
  phi = (kinParams.waveSpeed*t-s)/(kinParams.waveSpeed/kinParams.beatFreq);
  
  % Calculate phase
  L = phi - floor(phi);
  
  % Calculate inflection point age
  N = 2*phi - floor(2*phi);
  N = s/kinParams.waveSpeed + 0.5*N/kinParams.beatFreq;
  
  % Which have left/right curvature
  left = L < 0.5;
  right = ~left;
  
  % Calculate curvature function
  kL = 0.5*kinParams.kLeft*(1-cos(2*pi*N/kinParams.mLeft));
  kR = 0.5*kinParams.kRight*(1-cos(2*pi*N/kinParams.mRight));

  % Make curvature
  K = -kL.*left + kR.*right;
  
  % Calculate tail angle
  theta = kinParams.thetaAmp*cos(2*pi*kinParams.beatFreq*t) ...
    + kinParams.thetaOff;
  
  dtheta = ds.*K;
  theta = theta + cumsum(dtheta,2);
  
  dx = ds.*cos(theta);
  x = cumsum(dx,2);
  
  dy = ds.*sin(theta);
  y = cumsum(dy,2);

  % Filter out numerical junk
  [b,a] = butter(4, 1/(2*upsamp));
  
  for i=1:size(x,2)
    x(:,i) = filtfilt(b,a,x(:,i));
    y(:,i) = filtfilt(b,a,y(:,i));
  end
  
  % Downsample everything
  for i=1:size(x,2)
    s2(:,i) = downsample(s(:,i),upsamp);
    t2(:,i) = downsample(t(:,i),upsamp);
    x2(:,i) = downsample(x(:,i),upsamp);
    y2(:,i) = downsample(y(:,i),upsamp);
  end
  
  % Find points at center of plate
  s = (s2(:,1:end-1)+s2(:,2:end))/2;
  t = (t2(:,1:end-1)+t2(:,2:end))/2;
  x = -(x2(:,1:end-1)+x2(:,2:end))/2;
  y = (y2(:,1:end-1)+y2(:,2:end))/2;
  
  % Reduce s and t to row vector
  s = s(1,:);
  t = t(:,1)';
  
   
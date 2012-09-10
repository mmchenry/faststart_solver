function I = inertiaTensor(trunkParams,ocellParams,finParams,meatParams,tailX,tailY,Cm)
% This function calculates the center of mass,
% for the given morphological parameters
%
% Input data structure
%  trunkParams
%    -> s - row vector
%    -> right - row vector
%    -> left - row vector
%    -> dorsal - row vector
%    -> ventral - row vector
%    -> density - scalar
%
%  ocellParams
%    -> antPost - scalar
%    -> leftRight - scalar 
%    -> dorsoVent - scalar
%    -> radius - scalar
%    -> density - scalar
%
%  finParams
%    -> s - row vector
%    -> height - row vector
%    -> depth - row vector
%    -> width - row vector
%    -> density - scalar
%
%  meatParams
%    -> s - row vector
%    -> radius - scalar
%    -> density - scalar
%
%  tailX, tailY - row vector (same length as finParams.s &
%     meatParams.s), all start from x = -inf
%
%  Cm
%    -> x - scalar
%    -> y - scalar 
%    -> z - scalar
%
% Output data structure
%  I (1:3,1:3)
%  

  % Check assumptions
  % -----------------------------------------------------
  if (max(abs(trunkParams.right-trunkParams.left)) > 100*eps)
    error('trunk left and right must be equal')
  end

  % Calculate Trunk Contribution
  % -----------------------------------------------------
  a = (trunkParams.dorsal+trunkParams.ventral)/2;
  b = (trunkParams.right+trunkParams.left)/2;
  w = mean(diff(trunkParams.s));
  x0 = trunkParams.s-Cm.x;
  y0 = -Cm.y*ones(1,length(trunkParams.s));
  z0 = trunkParams.dorsal-a-Cm.z;

  ITrunk(1,1) = trunkParams.density*trapz(trunkParams.s,...
      a.*b.*((pi/4)*(a.^2+b.^2)+pi*(z0.^2+y0.^2)));
  ITrunk(1,2) = trunkParams.density*trapz(trunkParams.s,...
      a.*b.*(-pi*x0.*y0));
  ITrunk(1,3) = trunkParams.density*trapz(trunkParams.s,...
      a.*b.*(-pi*x0.*z0));

  ITrunk(2,1) = ITrunk(1,2);
  ITrunk(2,2) = trunkParams.density*trapz(trunkParams.s,...
      a.*b.*((pi/12)*(3*b.^2+w^2)+pi*(x0.^2+z0.^2)));
  ITrunk(2,3) = trunkParams.density*trapz(trunkParams.s,...
      a.*b.*(-pi*y0.*z0));

  ITrunk(3,1) = ITrunk(1,3);
  ITrunk(3,2) = ITrunk(2,3);
  ITrunk(3,3) = trunkParams.density*trapz(trunkParams.s,...
      a.*b.*((pi/12)*(3*a.^2+w^2)+pi*(x0.^2+y0.^2)));
  clear a b w x0 z0 y0;

  % Calculate Dorsal Tail Fin Contribution
  % -----------------------------------------------------
  w = finParams.width;
  h = finParams.height-meatParams.radius;
  L = mean(diff(finParams.s));
  x0 = tailX-Cm.x;
  y0 = tailY-Cm.y;
  z0 = (finParams.height-meatParams.radius)/2+meatParams.radius-Cm.z;

  IDorsalTailFin(1,1) = finParams.density*trapz(finParams.s,...
      h.*w.*((1/12)*(h.^2+w.^2)+y0.^2+z0.^2));
  IDorsalTailFin(1,2) = finParams.density*trapz(finParams.s,...
      h.*w.*(-x0.*y0));
  IDorsalTailFin(1,3) = finParams.density*trapz(finParams.s,...
      h.*w.*(-x0.*z0));

  IDorsalTailFin(2,1) = IDorsalTailFin(1,2);
  IDorsalTailFin(2,2) = finParams.density*trapz(finParams.s,...
      h.*w.*((1/12)*(L.^2+h.^2)+x0.^2+z0.^2));
  IDorsalTailFin(2,3) = finParams.density*trapz(finParams.s,...
      h.*w.*(-y0.*z0));

  IDorsalTailFin(3,1) = IDorsalTailFin(1,3);
  IDorsalTailFin(3,2) = IDorsalTailFin(2,3);
  IDorsalTailFin(3,3) = finParams.density*trapz(finParams.s,...
      h.*w.*((1/12)*(L.^2+w.^2)+x0.^2+y0.^2));
  clear w h L x0 y0 z0;

  % Calculate Ventral Tail Fin Contribution
  % -----------------------------------------------------
  w = finParams.width;
  h = finParams.depth-meatParams.radius;
  L = mean(diff(finParams.s));
  x0 = tailX-Cm.x;
  y0 = tailY-Cm.y;
  z0 = -((finParams.depth-meatParams.radius)/2+meatParams.radius)-Cm.z;

  IVentralTailFin(1,1) = finParams.density*trapz(finParams.s,...
      h.*w.*((1/12)*(h.^2+w.^2)+y0.^2+z0.^2));
  IVentralTailFin(1,2) = finParams.density*trapz(finParams.s,...
      h.*w.*(-x0.*y0));
  IVentralTailFin(1,3) = finParams.density*trapz(finParams.s,...
      h.*w.*(-x0.*z0));


  IVentralTailFin(2,1) = IVentralTailFin(1,2);
  IVentralTailFin(2,2) = finParams.density*trapz(finParams.s,...
      h.*w.*((1/12)*(L.^2+h.^2)+x0.^2+z0.^2));
  IVentralTailFin(2,3) = finParams.density*trapz(finParams.s,...
      h.*w.*(-y0.*z0));

  IVentralTailFin(3,1) = IVentralTailFin(1,3);
  IVentralTailFin(3,2) = IVentralTailFin(2,3);
  IVentralTailFin(3,3) = finParams.density*trapz(finParams.s,...
      h.*w.*((1/12)*(L.^2+w.^2)+x0.^2+y0.^2));
  clear w h L x0 y0 z0;

  % Calculate Tail Meat Contribution
  % -----------------------------------------------------
  r = meatParams.radius;
  L = mean(diff(meatParams.s));
  x0 = tailX-Cm.x;
  y0 = tailY-Cm.y;
  z0 = -Cm.z;

  ITailMeat(1,1) = pi*meatParams.density*trapz(meatParams.s,...
      r.^2.*((1/2)*(r.^2)+y0.^2+z0.^2));
  ITailMeat(1,2) = pi*meatParams.density*trapz(meatParams.s,...
      r.^2.*(-x0.*y0));
  ITailMeat(1,3) = pi*meatParams.density*trapz(meatParams.s,...
      r.^2.*(-x0.*z0));

  ITailMeat(2,1) = ITailMeat(1,2);
  ITailMeat(2,2) = pi*meatParams.density*trapz(meatParams.s,...
      r.^2.*((1/12)*(L.^2+3*(r.^2))+x0.^2+z0.^2));
  ITailMeat(2,3) = pi*meatParams.density*trapz(meatParams.s,...
      r.^2.*(-y0.*z0));

  ITailMeat(3,1) = ITailMeat(1,3);
  ITailMeat(3,2) = ITailMeat(2,3);
  ITailMeat(3,3) = pi*meatParams.density*trapz(meatParams.s,...
      r.^2.*((1/12)*(L.^2+3*r.^2)+x0.^2+y0.^2));
  clear r L x0 y0 z0;

  % Calculate Ocellus Contribution
  % -----------------------------------------------------
  r = ocellParams.radius;
  x0 = ocellParams.antPost-Cm.x;
  y0 = ocellParams.leftRight-Cm.y;
  z0 = ocellParams.dorsoVent-Cm.z;

  IOcellus(1,1) = (4*pi/3)*(ocellParams.density-trunkParams.density)*r^3*...
      ((2/5)*r^2+y0^2+z0^2);
  IOcellus(1,2) = (4*pi/3)*(ocellParams.density-trunkParams.density)*r^3*(-x0*y0);
  IOcellus(1,3) = (4*pi/3)*(ocellParams.density-trunkParams.density)*r^3*(-x0*z0);

  IOcellus(2,1) = IOcellus(1,2);
  IOcellus(2,2) = (4*pi/3)*(ocellParams.density-trunkParams.density)...
      *r^3*((2/5)*r^2+x0^2+z0^2);
  IOcellus(2,3) = (4*pi/3)*(ocellParams.density-trunkParams.density)*r^3*(-y0*z0);

  IOcellus(3,1) = IOcellus(1,3);
  IOcellus(3,2) = IOcellus(2,3);
  IOcellus(3,3) = (4*pi/3)*(ocellParams.density-trunkParams.density)...
      *r^3*((2/5)*r^2+x0^2+y0^2);
  clear r x0 y0 z0;

  % Big Addition
  % -----------------------------------------------------
  I = ITrunk + IDorsalTailFin + IVentralTailFin + ITailMeat + IOcellus;


function [M,x,y,z] = bodyMass(trunkParams,ocellParams,finParams,meatParams,tailX,tailY)
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
%    -> s - row vector (same length as tailX & tailY)
%    -> height - row vector
%    -> depth - row vector
%    -> width - row vector
%    -> density - scalar
%
%  meatParams
%    -> s - row vector (same length as tailX & tailY)
%    -> radius - row vector
%    -> density -scalar
%
%  tailX, tailY - row vector (same length as meatParams.s &
%     finParams.s), all start from x = -inf
  
%
% Output data structure
%  M - mass of the body
%  x,y,z - position of the center of mass

  
  % Check assumptions
  % -----------------------------------------------------
  if (max(abs(trunkParams.right-trunkParams.left)) > 100*eps)
    error('trunk left and right must be equal')
  end

  % Calculation of total mass
  % -----------------------------------------------------
  dorsalFinMass = finParams.density*trapz(finParams.s,...
      finParams.width.*(finParams.height-meatParams.radius));

  ventralFinMass = finParams.density*trapz(finParams.s,...
      finParams.width.*(finParams.depth-meatParams.radius));

  tailMeatMass = meatParams.density*trapz(meatParams.s,...
      pi*meatParams.radius.^2);

  trunkMass = trunkParams.density*trapz(trunkParams.s,...
      0.5*pi*trunkParams.right.*(trunkParams.dorsal+trunkParams.ventral))...
      -trunkParams.density*(4*pi/3)*ocellParams.radius^3;

  ocellusMass = ocellParams.density*(4*pi/3)*ocellParams.radius^3;

  M = dorsalFinMass + ventralFinMass + tailMeatMass + trunkMass + ocellusMass;

  % Calculation of the x position of the center of mass
  % -----------------------------------------------------  
  dorsalFinX = finParams.density*trapz(finParams.s,...
      finParams.width.*(finParams.height-meatParams.radius).*tailX);

  ventralFinX = finParams.density*trapz(finParams.s,...
      finParams.width.*(finParams.depth-meatParams.radius).*tailX);

  tailMeatX = meatParams.density*trapz(meatParams.s,...
      pi*meatParams.radius.^2.*tailX);

  trunkX = trunkParams.density*trapz(trunkParams.s,...
      0.5*pi*trunkParams.right.*(trunkParams.dorsal+trunkParams.ventral).*trunkParams.s)...
      -trunkParams.density*(4*pi/3)*ocellParams.radius^3*ocellParams.antPost;

  ocellusX = ocellParams.density*(4*pi/3)*ocellParams.radius^3*ocellParams.antPost;

  x = (dorsalFinX + ventralFinX + tailMeatX + trunkX + ocellusX)/M;

  % Calculation of the y position of the center of mass
  % -----------------------------------------------------
  dorsalFinY = finParams.density*trapz(finParams.s,...
      [finParams.width.*(finParams.height-meatParams.radius).*tailY]);
  
  ventralFinY= finParams.density*trapz(finParams.s,...
      [finParams.width.*(finParams.depth-meatParams.radius).*tailY]);

  tailMeatY = meatParams.density*trapz(meatParams.s,...
      pi*meatParams.radius.^2.*tailY);

  % By assumption of frontal symmetry
  trunkY = -trunkParams.density*(4*pi/3)*ocellParams.radius^3*ocellParams.leftRight;

  ocellusY = ocellParams.density*(4*pi/3)*ocellParams.radius^3*ocellParams.leftRight;

  y = (dorsalFinY + ventralFinY + tailMeatY + trunkY + ocellusY)/M;

  % Calculation of the z position of the center of mass
  % -----------------------------------------------------
  dorsalFinZ = finParams.density*trapz(finParams.s,...
      finParams.width.*(finParams.height-meatParams.radius).*...
      ((finParams.height-meatParams.radius)/2+meatParams.radius));
  
  ventralFinZ = -finParams.density*trapz(finParams.s,...
      finParams.width.*(finParams.depth-meatParams.radius).*...
      0.5.*(finParams.depth+meatParams.radius));

  tailMeatZ = 0;      % By assumption of circular cross-section

  trunkZ = trunkParams.density*trapz(trunkParams.s,...
      pi*trunkParams.right.*0.5.*(trunkParams.dorsal+trunkParams.ventral).*...
      0.5.*(trunkParams.dorsal-trunkParams.ventral)) ...
      -trunkParams.density*(4*pi/3)*ocellParams.radius^3*ocellParams.dorsoVent;

  ocellusZ = ocellParams.density*(4*pi/3)*ocellParams.radius^3*ocellParams.dorsoVent;

  z = (dorsalFinZ + ventralFinZ + tailMeatZ + trunkZ + ocellusZ)/M;

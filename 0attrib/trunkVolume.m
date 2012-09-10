function [V,x,y,z] = trunkVolume(trunkParams)
% This function calculates the trunk center of volume,
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
% Output data structure
%  V - volume of the trunk
%  x,y,z - position of the trunk center of volume

  % Check assumptions
  % -----------------------------------------------------
  if (max(abs(trunkParams.right-trunkParams.left)) > 100*eps)
    error('trunk left and right must be equal')
  end

  % Calculation of total volume
  % -----------------------------------------------------
  trunkVol = trapz(trunkParams.s, pi*trunkParams.right...
      .*((trunkParams.dorsal+trunkParams.ventral)/2));
  V = trunkVol;

  % Calculation of the x position of the center of volume
  % -----------------------------------------------------
  trunkX = trapz(trunkParams.s,pi*trunkParams.right...
      .*((trunkParams.dorsal+trunkParams.ventral)/2).*trunkParams.s);
  x = trunkX/V;

  % Calculation of the y position of the center of volume
  % -----------------------------------------------------
  trunkY = 0;    % By assumption of frontal symmetry
  y = trunkY/V;

  % Calculation of the z position of the center of volume
  % -----------------------------------------------------
  trunkZ = trapz(trunkParams.s,pi*trunkParams.right...
      .*((trunkParams.dorsal+trunkParams.ventral)/2)...
      .*((trunkParams.dorsal-trunkParams.ventral)/2));
  z = trunkZ/V;


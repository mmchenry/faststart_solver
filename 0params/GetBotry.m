function [trunkParams,finParams,meatParams,ocellParams,kinParams] = getBotry(simParams,ROOT)
% This function returns all of the parameters that describe the 
% morphology and kinematics of Botrylloides spp.
%
% Input data structure
%  simParams
%    -> numTailSegs - scalar
%    -> numTrunkSegs - scalar
%    -> numTailBeats - scalar
%    -> sampleRate - scalar
%  
% Output data structure
%  trunkParams
%    -> s - row vector
%    -> right - row vector
%    -> left - row vector
%    -> dorsal - row vector
%    -> ventral - row vector
%    -> density - scalar
%
%  ocellParams
%    -> direction - column vector (3x1)
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
%    -> radius - row vector
%    -> density - scalar
%
%  kinParams
%    -> kLeft - scalar
%    -> kRight - scalar
%    -> mLeft - scalar
%    -> mRight - scalar
%    -> thetaAmp - scalar
%    -> thetaOff - scalar
%    -> waveSpeed - scalar
%    -> beatFreq - scalar
%
% Units: 10^-4 m, 10^-2 s, 10^-4 g
  
  % Set file locations
  % -----------------------------------------------------    
  trunk_morpho_file =   [ROOT filesep '0params' filesep '0botry' filesep 'trunk_morpho.mat'];
  ocell_morpho_file =   [ROOT filesep '0params' filesep '0botry' filesep 'ocell_morpho.mat'];
  fin_morpho_file =     [ROOT filesep '0params' filesep '0botry' filesep 'fin_morpho.mat'];
  meat_morpho_file =    [ROOT filesep '0params' filesep '0botry' filesep 'meat_morpho.mat'];
  kine_morpho_file =    [ROOT filesep '0params' filesep '0botry' filesep 'kine_morpho.mat'];

  % Put together trunk data
  % -----------------------------------------------------    
  load (trunk_morpho_file);
  trunkParams.s = linspace(min(trunkData.s),max(trunkData.s),simParams.numTrunkSegs);
  trunkParams.right = csapi(trunkData.s,trunkData.right,trunkParams.s);
  trunkParams.left = csapi(trunkData.s,trunkData.left,trunkParams.s);
  trunkParams.dorsal = csapi(trunkData.s,trunkData.dorsal,trunkParams.s);
  trunkParams.ventral = csapi(trunkData.s,trunkData.ventral,trunkParams.s);
  trunkParams.density = trunkData.density;
  clear trunkData

  % Put together ocellus data
  % -----------------------------------------------------  
  load (ocell_morpho_file);
  ocellParams.direction = ocellData.direction;
  ocellParams.antPost = ocellData.antPost;
  ocellParams.leftRight = ocellData.leftRight;
  ocellParams.dorsoVent = ocellData.dorsoVent;
  ocellParams.radius = ocellData.radius;
  ocellParams.density = ocellData.density;
  clear ocellData;

  % Put together tail fin data
  % -----------------------------------------------------  
  load (fin_morpho_file); 
  finParams.s = linspace(min(finData.s),max(finData.s),simParams.numTailSegs);
  finParams.height = csapi(finData.s,finData.height,finParams.s);
  finParams.depth = csapi(finData.s,finData.depth,finParams.s);
  finParams.width = csapi(finData.s,finData.width,finParams.s);
  finParams.density = finData.density;
  clear finData;
 
  % Put together tail meat data
  % -----------------------------------------------------  
  load (meat_morpho_file); 
  meatParams.s = linspace(min(meatData.s),max(meatData.s),simParams.numTailSegs);
  meatParams.radius = csapi(meatData.s,meatData.radius,meatParams.s);
  meatParams.density = meatData.density;
  clear meatData;
 
  % Put together kinematic data
  % -----------------------------------------------------
  load (kine_morpho_file);
  kinParams.kLeft = kineData.kLeft;
  kinParams.kRight = kineData.kRight;
  kinParams.mLeft = kineData.mLeft;
  kinParams.mRight = kineData.mRight;
  kinParams.thetaAmp = kineData.thetaAmp;
  kinParams.thetaOff = kineData.thetaOff;
  kinParams.waveSpeed = kineData.waveSpeed;
  kinParams.beatFreq = kineData.beatFreq;
  clear kineData;
  
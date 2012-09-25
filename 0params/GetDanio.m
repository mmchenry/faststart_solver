function [finParams,meatParams,ocellParams,kinParams] = GetBotry(simParams,ROOT)
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
%  kinParams.stg1
%    -> dur - scalar
%    -> kMax - scalar
%
%  kinParams.stg2
%    -> kLeft - scalar
%    -> kRight - scalar
%    -> mLeft - scalar
%    -> mRight - scalar
%    -> waveSpeed - scalar
%    -> beatPeriod - scalar
%
%  kinParams.und
%    -> kLeft - scalar
%    -> kRight - scalar
%    -> mLeft - scalar
%    -> mRight - scalar
%    -> waveSpeed - scalar
%    -> beatPeriod - scalar
%
% Units: 10^-4 m, 10^-2 s, 10^-4 g

%% Set file locations
ocell_morpho_file =   [ROOT filesep '0params' filesep '0danio' filesep 'ocell_morpho.mat'];
fin_morpho_file =     [ROOT filesep '0params' filesep '0danio' filesep 'fin_morpho.mat'];
meat_morpho_file =    [ROOT filesep '0params' filesep '0danio' filesep 'meat_morpho.mat'];
kine_morpho_file =    [ROOT filesep '0params' filesep '0danio' filesep 'kine_morpho.mat'];


%% Put together ocellus data

load (ocell_morpho_file);
ocellParams.direction = ocellData.direction;
ocellParams.antPost = ocellData.antPost;
ocellParams.leftRight = ocellData.leftRight;
ocellParams.dorsoVent = ocellData.dorsoVent;
ocellParams.radius = ocellData.radius;
ocellParams.density = ocellData.density;

clear ocellData;


%% Put together tail fin data

load (fin_morpho_file);
finParams.s = linspace(min(finData.s),max(finData.s),simParams.numTailSegs);
finParams.height = csapi(finData.s,finData.height,finParams.s);
finParams.depth = csapi(finData.s,finData.depth,finParams.s);
finParams.width = csapi(finData.s,finData.width,finParams.s);
finParams.density = finData.density;

clear finData;


%% Put together tail meat data

load (meat_morpho_file);
meatParams.s = linspace(min(meatData.s),max(meatData.s),simParams.numTailSegs);
meatParams.radius = csapi(meatData.s,meatData.radius,meatParams.s);
meatParams.density = meatData.density;

clear meatData;


%% Put together kinematic data

load (kine_morpho_file);

% Body position where bending begins
kinParams.s_startBend = kineData.s_startBend;

% Parameters for stage 1 of fast start -------------------------------

% Stage 1 duration
kinParams.stg1.dur = kineData.stg1.dur;

% Wave speed
kinParams.stg1.waveSpeed = kineData.stg1.waveSpeed;

% Max curvature
kinParams.stg1.kMax = kineData.stg1.kMax;


% Parameters for stage 2 of fast start -------------------------------

% Wave speed & duration (half beat)
kinParams.stg2.dur       = kineData.stg2.beatPeriod;
kinParams.stg2.waveSpeed = kineData.stg2.waveSpeed;

% Max curvature on left and right sides
kinParams.stg2.kLeft     = kineData.stg2.kLeft;
kinParams.stg2.kRight    = kineData.stg2.kRight;


% Parameters for undulatory swimming ---------------------------------

% Wave speed & beat period (full beat cycle)
kinParams.und.beatPeriod = kineData.und.beatPeriod;
kinParams.und.waveSpeed  = kineData.und.waveSpeed;

% Max curvature on left and right sides
kinParams.und.kLeft     = kineData.und.kLeft;
kinParams.und.kRight    = kineData.und.kRight;

clear kineData;

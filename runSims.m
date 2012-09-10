function runSims()
% Function collects parameters from file
% and then runs through a parameter space

% Set root
ROOT = '/Volumes/Docs/Dropbox/matlab/faststart_model/faststart_solver';

% Set output
outfn = [ROOT filesep 'data.mat'];

% Source it up
addpath([ROOT filesep '0attrib']);
addpath([ROOT filesep '0forces']);
addpath([ROOT filesep '0params']);
addpath([ROOT filesep '0trans']);

% Prepare standard
p = GetParams('botry',ROOT);
%p = GetParams('default');

% Add behavioral parameters
p.behaveModel = 1;
p.behaveGain = 0;
p.behaveAmp = pi/4;
p.behaveThresh = .5;
p.behaveLag = .5;
p.behaveOffset = 0;

% Clear global brightness
global gBright gBrightT
gBright = [];
gBrightT = [];

% Solve
data = solver(p);
data.p = p;

% Save
save(outfn,'data')

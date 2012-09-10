function runSims()
% Function collects parameters from file
% and then runs through a parameter space

  % Set root
  %ROOT = '/home/jims/ascidian/model/solver';
  ROOT = '/Users/mmchenry/Dropbox/Matlab/faststart_model/solver';

  % Set output
  outfn = [ROOT filesep 'data.mat'];
  
  % Source it up
  addpath([ROOT filesep '0attrib']);
  addpath([ROOT filesep '0forces']);
  addpath([ROOT filesep '0params']);
  addpath([ROOT filesep '0trans']);

  % Prepare standard
  p = GetParams('default');

  % Modify behavior
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

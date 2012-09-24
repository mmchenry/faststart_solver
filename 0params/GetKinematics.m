function [t,s,x,y] = GetKinematics(kinParams,simParams,tailLength,timeStop)
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
%  timeStop   - scalar
%
% Output data structure
%  t,s,x,y are matrices with the arc length
%  in rows and the time in columns
%

% Set oversample factor
upsamp = 10;

% Generate arc length, time grid
[s,t] = meshgrid(linspace(0,tailLength,simParams.numTailSegs+1),...
                 linspace(0,timeStop,upsamp*simParams.sampleRate));
             
% Positon where body starts to bend            
s_start  = kinParams.s_startBend;


% TODO: Add rigid trunk region
% TODO: Adjust time to correct values 
% TODO: delete unused kinematic parameters & correct time values
% TODO: make change in period from stage 2 to 3 more realistic 
% (check Muller)


%% Stage 1    
% The body curls into a "C" toward the left
     
% Parameters
wv_spd   = kinParams.stg1.waveSpeed;
b_period = kinParams.stg1.dur;
kLeft  = -kinParams.stg1.kMax;
kRight = 0;
stg1_end = (b_period + s(end)/wv_spd);

% Index of time for duration of stage 1
t_start = t(1);
t_end   = t(find(...
            t(:,1)<=stg1_end,...
                     1,'last'),1);
t_idx   = (t>=t_start) & (t<=t_end);

% Wave position
phi1 = t_idx .* (wv_spd*t-(s-0*s_start))/(wv_spd*b_period);

% Curvature
K1 = -((kLeft-kRight).*(0.5.*cos(pi.*phi1./b_period)) - ((kLeft+kRight)./2));

% Define position-specifc start and end of bending wave
wv_start = s/wv_spd;
wv_end = b_period+s/wv_spd;

% Set curvature constant outside of wave propagation
K1(t>wv_end)   = kLeft;
K1(t<wv_start) = kRight;

% Zero out values outside of stage 1
K1 = t_idx .* K1;

% Set values before delay to stage 1 curvature
%kLast = repmat(K1(find(t_idx(:,1)==1,1,'last'),:),size(t,1),1);

clear wv_spd b_period kLeft kRight


%% Stage 2 
% The body unfurls.  It ends when the anterior-most point reaches the
% maximum right curvature.

% Parameters
wv_spd      = kinParams.stg2.waveSpeed;
b_period    = kinParams.stg2.beatPeriod;
kLeft       = -kinParams.stg1.kMax;
kRight      = kinParams.stg2.kRight;

stg2_end    = stg1_end + b_period;

% Index of time for duration of stage 2
t_start = t(find(...
            t(:,1)>stg1_end,...
                     1,'first'),1);
t_end   = t(find(...
            t(:,1)<=stg2_end,...
                     1,'last'),1);
t_idx   = (t>=t_start) & (t<=t_end);


%b_period = repmat(linspace(b_pd_min,b_pd_max,size(t,2)),size(t,1),1);

% Calculate wave position
phi2 = (wv_spd*(t-t_start)-(s-0*s_start))./(wv_spd*b_period);

% Calculate curvature
K2 = ((kLeft-kRight).*(0.5.*cos(pi.*phi2)) + ((kLeft+kRight)./2));

% Define delay for the start of bending (varies with s)
delay2 = s/wv_spd;

% Set values before delay to stage 1 curvature
%kLast = repmat(K1(end,:),size(t,1),1);
K2((t-t_start)<delay2) = kLeft;

%K2((t-t_start)<delay) = kLast((t-t_start)<delay);

% Zero out values outside of stage 2
K2 = t_idx .* K2;

% Calculate position at end of this phase
i_last = find(t_idx(:,1)==1,1,'last');
end_phs = (max(delay2(:))-delay2)./b_period;

clear wv_spd b_period kLeft kRight t_idx t_start t_end


%% Stage 3 
% Transition into regular undulation

% Define short names
wv_spd      = kinParams.stg2.waveSpeed;
b_pd_min    = kinParams.stg2.beatPeriod;
b_pd_max    = kinParams.und.beatPeriod;
kLeft       = -kinParams.stg2.kLeft;
kRight      = kinParams.stg2.kRight;
stg3_end    = stg2_end + b_pd_min;

b_period = repmat(linspace(b_pd_min,b_pd_max,size(t,2)),size(t,1),1);

% Index of time for duration of stage 3
t_start = t(find(...
            t(:,1)>stg2_end,...
                     1,'first'),1);
t_end   = t(find(...
            t(:,1)<=stg3_end,...
                     1,'last'),1);
t_idx   = (t>=t_start) & (t<=t_end);

% Calculate wave position
phi3 = (wv_spd*(t-t_start)-(s-0*s_start))./(wv_spd*b_period);

% Stage 3 curvature
K3 = ((kRight-kLeft).*(0.5.*cos(pi.*phi3)) + ((-kLeft-kRight)./2));

% Zero out values outside of stage 3
K3 = t_idx .* K3;

clear wv_spd b_period kLeft kRight t_idx t_start t_end


%% Stage 4
% Transition into regular undulation

% Define short names
wv_spd      = kinParams.und.waveSpeed;
b_period    = kinParams.und.beatPeriod;
kLeft       = -kinParams.und.kLeft;
kRight      = kinParams.und.kRight;

stg4_end    = stg3_end + 2*b_period.*simParams.numTailBeats;

% Index of time for duration of stage 2
t_start = t(find(...
            t(:,1)>stg3_end,...
                     1,'first'),1);
t_end   = t(find(...
            t(:,1)<=stg4_end,...
                     1,'last'),1);
t_idx   = (t>=t_start) & (t<=t_end);

% Calculate wave position
phi4 = (wv_spd*(t-t_start)-(s-0*s_start))...
    /(wv_spd*b_period);

% Stage 4 curvature
K4 = ((kRight-kLeft).*(0.5.*cos(pi.*phi3)) + ((-kLeft-kRight)./2));

% Zero out values outside of stage 3
K4 = t_idx .* K4;



%% Combine all stages
K_tot = K1 + K2 + K3 + K4;

% Define delay for the start of bending
delay = s/wv_spd;

K2((t-t_start)<delay) = kLeft;

K2 = t_idx .* K2;

%k = k .* (s>s_start);




for i = 1:t_tmp
    
   
    
    
end


% Calculate wave position
phi = (wv_spd*t2-s)/(wv_spd*b_period);

%t2 = t(t_idx);

% Calculate phase
L = phi - floor(phi);

% Calculate inflection point age
N = 2*phi - floor(2*phi);
N = s/wv_spd + 0.5*N*b_period;

% Which have left/right curvature
left = L < 0.5;
right = ~left;

% Calculate curvature function
kL = 0.5*kLeft*(1-cos(2*pi*N/mLeft));
kR = 0.5*kRight*(1-cos(2*pi*N/mRight));

% Make curvature
K2 = -kL.*left + kR.*right;


% Set stage 2 curvature equal to values at end of stage 1 
K2 = K(find(t(:,1)<=kinParams.stg1.dur,1,'last'),:);
K2 = repmat(K2,max(sum(t_idx(:,1))),1);


% Clear short names
clear wv_spd b_period kLeft kRight right left s_start t_start t_end

% % Curvature during stage 2
% K_stg2 = K_stage2(s,t,kinParams.s_startBend,kinParams.stg1.dur, ...
%                   kinParams.stg2.waveSpeed,kinParams.und.kLeft, ...
%                   kinParams.und.kRight,kinParams.und.mLeft,...
%                   kinParams.und.mRight,K_stg1(end,:));
% 
% % Curvature during undulation
% K = K_undulation(s,t,kinParams.und.waveSpeed,kinParams.und.beatPeriod,...
%             kinParams.und.kLeft,kinParams.und.kRight,kinParams.und.mLeft,...
%             kinParams.und.mRight);


% Trim region of zero curvature (i.e. within the cranium)
K = K .* (s>=kinParams.s_startBend);


animateK(s,t,K)
% Calculate ds
ds = diff(s,1,2);
ds = [ds(:,1) ds];

% Calculate tail angle from curvature
dtheta = ds.*K;
theta = cumsum(dtheta,2);

% Calculate x and y positions from angle
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




function K_stg2 = K_stage2(s,t,s_start,stg1_dur,waveSpeed,kLeft,kRight,mLeft,...
                           mRight,k_last)
                       
% TODO: Fix Stage 2 (needs to transition from stage 1 to undulatory mode, 
% over the duration necessary to have a wave travel down the body.                      
                       
                       
% Returns matrix of curvature values for fast start stage 2

% Stage 2 duration (time to pass 1 wave)
stg2_dur = (max(s(:))-s_start) / waveSpeed;

% Time values for stage 2 
t_curr = t((t(:,1)>=stg1_dur) & (t(:,1) < stg2_dur),1);

K_stg2 = k_last;

% Loop thru time 
for i = 2:length(t_curr)
    
    
    %TODO: Fix below
    K_stg2(i,:) = K_stg2(i-1,:) + 2;
    
    
end

% Trim values for first time
K_stg2 = K_stg2(2:end,:);

% Calculate wave position
phi = (waveSpeed*t-s)/(waveSpeed*beatPeriod);

% Calculate phase
L = phi - floor(phi);

% Calculate inflection point age
N = 2*phi - floor(2*phi);
N = s/waveSpeed + 0.5*N*beatPeriod;

% Which have left/right curvature
left = L < 0.5;
right = ~left;

% Calculate curvature function
kL = 0.5*kLeft*(1-cos(2*pi*N/mLeft));
kR = 0.5*kRight*(1-cos(2*pi*N/mRight));

% Make curvature
K = -kL.*left + kR.*right;


function K = K_undulation(s,t,waveSpeed,beatPeriod,kLeft,kRight,mLeft,mRight)
% Returns matrix of curvature values for undulatory swimming

% Calculate wave position
phi = (waveSpeed*t-s)/(waveSpeed*beatPeriod);

% Calculate phase
L = phi - floor(phi);

% Calculate inflection point age
N = 2*phi - floor(2*phi);
N = s/waveSpeed + 0.5*N*beatPeriod;

% Which have left/right curvature
left = L < 0.5;
right = ~left;

% Calculate curvature function
kL = 0.5*kLeft*(1-cos(2*pi*N/mLeft));
kR = 0.5*kRight*(1-cos(2*pi*N/mRight));

% Make curvature
K = -kL.*left + kR.*right;

function animateK(s,t,K)

% Don't animate times with no curvature
idx = find(~(K(:,1)==0));

% Calculate ds
ds = diff(s,1,2);
ds = [ds(:,1) ds];

% Calculate tail angle from curvature
dtheta = ds.*K;
theta = cumsum(dtheta,2);

% Calculate x and y positions from angle
dx = ds.*cos(theta);
x = cumsum(dx,2);

dy = ds.*sin(theta);
y = cumsum(dy,2);


f = figure;
set(f,'DoubleBuffer','on');

for i = 1:length(idx)
    
    plot(x(idx(i),:),y(idx(i),:))
    axis equal
    pause(.001)
end


function plot_pcolor(s,t,M)

% Exclude zero values
idx = find(min(M,[],2)~=0);

% Plot and tweak
h = pcolor(s(idx,:),t(idx,:),M(idx,:));
xlabel('arclength')
ylabel('time')
set(h,'EdgeColor','none')
colorbar






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

%Indicies for bending region of body
i_start = find(s(1,:)>s_start);


%% Stage 1    
% The body curls into a "C" toward the left with a rapid wavespeed
     
% Parameters
wv_spd     = kinParams.stg1.waveSpeed;
b_period   = kinParams.stg1.dur;
kLeft      = -kinParams.stg1.kMax;
kRight     = 0;

% Index of time for duration of stage 1
t_start1 = (s-s_start)/wv_spd;
t_end1   = t_start1 + b_period;
t_idx1   = t<t_end1;

%t_start = t(1);
%t_end   = t(find(t(:,1)<=stg1_end,1,'last'),1);

% Wave position
phi1 = (wv_spd*t-(s-s_start))/(wv_spd*b_period);

% Curvature
K1 = -((kLeft-kRight).*(0.5.*cos(pi.*phi1./b_period)) - ((kLeft+kRight)./2));

% Define position-specifc start and end of bending wave
%wv_start = (s-s_start)/wv_spd;
%wv_end = b_period+(s-s_start)/wv_spd;

% Set to kRight for before wave propagation
%K1(t>b_period)   = kLeft;
K1(t<t_start1) = kRight;

% Make trunk/cranium rigid
K1(s<s_start) = 0;

% Zero out values outside of stage 1
K1 = t_idx1 .* K1;

% Clear out temporary variables
clear wv_spd b_period kLeft kRight wv_start wv_end t_idx1 t_start1


%% Stage 2 
% The body unfurls.  It ends when the anterior-most point reaches the
% maximum right curvature.

% Parameters
wv_spd      = kinParams.stg2.waveSpeed;
b_period    = kinParams.stg2.dur;
kLeft       = -kinParams.stg1.kMax;
kRight      = kinParams.stg2.kRight;
%stg2_end    = stg1_end + b_period;

% Define delay for the start of bending (varies with s)
delay2 = (s-s_start)/wv_spd;

% Index of time for duration of stage 2
t_start2 = min(t_end1(:,min(i_start))) + delay2;
t_end2   = t_start2 + b_period + delay2;
t_idx2   = (t>=t_end1) & (t<t_end2);

% Calculate wave position
phi2 = (wv_spd*(t-t_start2)-(s-s_start))./(wv_spd*b_period);

% Calculate curvature
K2 = ((kLeft-kRight).*(0.5.*cos(pi.*phi2)) + ((kLeft+kRight)./2));

% Set values before delay to stage 1 curvature
K2((t-t_start2)<delay2) = kLeft;

% Make trunk/cranium rigid
K2(s<s_start) = 0;

% Zero out values outside of stage 2
K2 = t_idx2 .* K2;

% Clear out temporary variables
clear wv_spd b_period kLeft kRight t_idx2 t_start2 delay2 t_start2
clear t_idx2


%% Stage 3 
% Transition into regular undulation

% Define short names
wv_spd      = kinParams.stg2.waveSpeed;
b_pd_min    = kinParams.stg2.dur;
b_pd_max    = kinParams.und.beatPeriod/2;
kLeft       = -kinParams.stg2.kLeft;
kRight      = kinParams.stg2.kRight;

%Define position-variable period
b_period = ones(size(t,1),size(t,2));
b_period(:,i_start) = ...
    repmat(linspace(b_pd_min,b_pd_max,length(i_start)),size(t,1),1);

% Define delay for the start of bending (varies with s)
delay3 = (s-s_start)/wv_spd;

% Index of time for duration of stage 2
t_start3 = min(t_end2(:,min(i_start))) + delay3;
t_end3   = t_start3 + b_period + delay3;
t_idx3   = (t>=t_end2) & (t<t_end3);

% Calculate wave position
phi3 = (wv_spd*(t-t_start3)-(s-s_start))./(wv_spd*b_period);

% Stage 3 curvature
K3 = ((kRight-kLeft).*(0.5.*cos(pi.*phi3)) + ((-kLeft-kRight)./2));

% Set values before delay to initial curvature
K3((t-t_start3)<delay3) = kRight;

% Make trunk/cranium rigid
K3(s<s_start) = 0;

% Zero out values outside of stage 3
K3 = t_idx3 .* K3;

% Clear out temporary variables
clear wv_spd b_period kLeft kRight t_idx3 t_start3 t_end2 phi3 b_pd_min b_pd_max


%% Stage 4 (undulation)
% Transition into regular undulation

% Define short names
wv_spd      = kinParams.und.waveSpeed;
b_period    = kinParams.und.beatPeriod;
kLeft       = -kinParams.und.kLeft;
kRight      = kinParams.und.kRight;

% Define delay for the start of bending (varies with s)
delay4 = (s-s_start)/wv_spd;

% Index of time for duration of stage 2
t_start4 = min(t_end3(:,min(i_start)))+ delay4;
t_end4   = min(t_start4(:)) + b_period.*simParams.numTailBeats;
t_idx4   = (t>=t_end3) & (t<t_end4);

% Calculate wave position
phi4 = (wv_spd*(t-t_start4)-(s-s_start))./(wv_spd*b_period/2);

% Stage 4 curvature
K4 = ((kLeft-kRight).*(0.5.*cos(pi.*phi4)) + ((kLeft+kRight)./2));

% Set values before delay to initial curvature
K4((t-t_start4)<delay4) = -kinParams.stg2.kLeft;

% Make trunk/cranium rigid
K4(s<s_start) = 0;

% Zero out values outside of stage 3
K4 = t_idx4 .* K4;

% Clear out temporary variables
clear wv_spd b_period kLeft kRight t_start t_end t_idx phi4


%% Combine all stages

% Add them together
K = K1 + K2 + K3 + K4;


%% Visualize curvature data
% Makes a pseudocolor plot of curvature, as a function of s and t

% Idicies for positions to examine
s1 = 55;
s2 = 122;
s3 = 190;

% Exclude zero values
idx = find(min(K,[],2)~=0);

% New window
figure

% Plot and tweak
subplot(4,1,[1:3])
h = pcolor(s(idx,:),t(idx,:),K(idx,:));
xlabel('arclength')
ylabel('time')
set(h,'EdgeColor','none')
colorbar
colormap cool

hold on
h1 = plot([s(1,s1) s(1,s1)],ylim,'k--');
h2 = plot([s(1,s2) s(1,s2)],ylim,'b--');
h3 = plot([s(1,s3) s(1,s3)],ylim,'r--');
set(h1,'LineWidth',1.5)
set(h2,'LineWidth',1.5)
set(h3,'LineWidth',1.5)

subplot(4,1,4)
h = plot(t(:,s1),K(:,s1),'k-',t(:,s2),K(:,s2),'b-',t(:,s3),K(:,s3),'r-');
xlabel('time')
ylabel('K')

clear h1 h2 h3 h s1 s2 s3 

% Animate midline
%figure
%animateK(s,t,K)


%% Post-processing

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
idx = find(~(K(:,end)==0));

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






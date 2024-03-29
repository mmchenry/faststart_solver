function kine_play

num_s = 20;

num_t = 50;

% Duration (s)
t_max = 20e-3;

% Body length (m)
s_max = 4e-3;

prd = 10e-3;

% Body position where bending starts (m)
pos_bend = 0.5e-3;

% Stage 1 duration (s)
stg1_dur = 10e-3;

% Stage 1 min_curvature (rad/m)
stg1_kappa_min = -2000;

% Stage 2 min_curvature (rad/m)
stg2_kappa_max = -1500;

% Stage 2 wavespeed (m/s)
stg2_wavespd = 0.8;

% Arclength (au)
s = linspace(0,s_max,num_s);

% Time (s)
t = linspace(0,t_max,num_t);

% Idicies
i_stg1 = t<=stg1_dur;
i_ant = s<=pos_bend;
i_post = s > pos_bend;

% Stage 1 values of over time   
stg1_kappa = repmat([stg1_kappa_min.*t(i_stg1)./stg1_dur]',1,length(s(i_post)));

% Stage 1 curvature values
kappa = [zeros(length(t(i_stg1)),length(s(i_ant))) stg1_kappa];
    
% Format position and time, compatible with kappa
[S,T] = meshgrid(s,t(i_stg1));


%TODO: Stage 2 kinematics

figure
surf(S,T,kappa)
view(2)
rotate3d
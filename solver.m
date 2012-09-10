function data = solver(p)
%  This function is the workhorse for the ascidian simulations.
%  It takes all of the parameters from master, and spits them
%  into a numerical integration function.
%
% Input data structure
%  p
%    -> timeStart - scalar
%    -> timeStop - scalar
%
%    -> light - column vector x N-time
%    -> gravity - column vector
%
%    -> fluidDensity - scalar
%    -> fluidKVisc - scalar
%
%    -> larvaTrunkVol - scalar
%    -> larvaTrunkCV - column vector
%
%    -> larvaTailT
%    -> larvaTailS
%    -> larvaTailRX
%    -> larvaTailRY
%    -> larvaTailUX
%    -> larvaTailUV
%    -> larvaTailAX
%    -> larvaTailAY
%
%    -> behaveModel - scalar
%    -> behaveGain - scalar
%    -> behaveAmp - scalar
%    -> behaveThresh - scalar
%    -> behaveLag - scalar
%
%    -> trunkParams - struct
%    -> ocellParams - struct
%    -> finParams- struct
%    -> meatParams - struct
%    -> kinParams - struct
%
% Output data structure
%    -> t
%    -> x
%    -> u
%    -> w
%    -> xi
%    -> theta
%    -> phi
%

% Setup initial conditions
% -----------------------------------------------------
xyz_initial     = [0 0 0];      % initial position in inertial frame
U_initial       = [0 5 0];      % initial velocity in inertial frame
w_initial       = [0 0 0];      % initial body rotational velocity
Euler_initial   = [0 pi/4 0];      % initial Euler angles  in xyz conv.

state0    = [xyz_initial U_initial w_initial Euler_initial];
timeSpan  = [p.timeStart p.timeStop];

% Setup ode options
% -----------------------------------------------------

options   = odeset('InitialStep', .005,...
    'RelTol',       1e-2,...
    'MaxStep',     .03,...
    'OutputFcn',    @DisplayData,...
    'Refine',       2,...
    'Vectorized',   'off');

% Run solver, Run !!!
% -----------------------------------------------------
[time state] = ode113(@EquationsOfMotion,timeSpan,state0,options,p);

% Save the data
% -----------------------------------------------------
data.t      = time';
data.x      = state(:,1:3)';
data.u      = state(:,4:6)';
data.w      = state(:,7:9)';
data.xi     = state(:,10)';
data.theta  = state(:,11)';
data.phi    = state(:,12)';

global gBright gBrightT
data.bright = gBright;





%------------------------------------------------------------------------------------
function f = EquationsOfMotion(t,stateMat,p)
%  Describes the equations of motion for the body of a tadpole in inertial coordinates.
%  state= [x y z Ux Uy Uz wx wy wz xi theta phi];

if (t > 80)
    disp('');
end

% Set state
state.t = t;
state.r = stateMat(1:3);
state.u = stateMat(4:6);
state.w = stateMat(7:9);
state.xi = stateMat(10);
state.theta = stateMat(11);
state.phi = stateMat(12);

% Set light intensity
state.light(1,1) = interp1(p.lightT,p.light(1,:),t);
state.light(2,1) = interp1(p.lightT,p.light(2,:),t);
state.light(3,1) = interp1(p.lightT,p.light(3,:),t);

state.dlight(1,1) = interp1(p.lightT,p.dlight(1,:),t);
state.dlight(2,1) = interp1(p.lightT,p.dlight(2,:),t);
state.dlight(3,1) = interp1(p.lightT,p.dlight(3,:),t);

state.light  = InertialToBody(state.light,state.xi,state.theta,state.phi);
state.dlight = InertialToBody(state.dlight,state.xi,state.theta,state.phi);

% Set percieved intensity
global gBright gBrightT
if (isempty(gBrightT) ||...
        min(gBrightT) > t-p.behaveLag ||...
        max(gBrightT) < t-p.behaveLag)
    state.bright = 0;
    state.dbright = 0;
else
    state.bright = interp1(gBrightT,gBright,t-p.behaveLag);
    state.dbright = interp1(gBrightT,[0 diff(gBright)],t-p.behaveLag);
end

% Set tail angle
switch p.behaveModel
    case 1
        state.tailAngle = p.behaveGain*state.bright+p.behaveOffset;
    case 2
        state.tailAngle = p.behaveGain*state.dbright+p.behaveOffset;
    case 3
        if (state.bright < -p.behaveThresh)
            state.tailAngle = - p.behaveAmp + p.behaveOffset;
        elseif (state.bright > p.behaveThresh)
            state.tailAngle = p.behaveAmp + p.behaveOffset;
        else
            state.tailAngle = p.behaveOffset;
        end
    case 4
        if (state.dbright < -p.behaveThresh)
            state.tailAngle = - p.behaveAmp + p.behaveOffset;
        elseif (state.dbright > p.behaveThresh)
            state.tailAngle = p.behaveAmp + p.behaveOffset;
        else
            state.tailAngle = p.behaveOffset;
        end
end

% Set tail points
nT = state.t - floor(state.t*p.kinParams.beatFreq)/p.kinParams.beatFreq;
TR(1,:) = interp2(p.larvaTailS,p.larvaTailT,...
                  p.larvaTailRX,p.larvaTailS(1,:),nT);
TR(2,:) = interp2(p.larvaTailS,p.larvaTailT,...
                  p.larvaTailRY,p.larvaTailS(1,:),nT);

rot = [cos(state.tailAngle)  sin(state.tailAngle);
       -sin(state.tailAngle) cos(state.tailAngle)];
TR = rot*TR;

state.tailrx = TR(1,:);
state.tailry = TR(2,:);


% plot(TR(1,:),TR(2,:))
% axis equal
% pause(.01)

clear TR rot nT

% Set tail velocity
nT = state.t - floor(state.t*p.kinParams.beatFreq)/p.kinParams.beatFreq;
state.tailux = interp2(p.larvaTailS,p.larvaTailT,...
    p.larvaTailUX,p.larvaTailS(1,:),nT);
state.tailuy = interp2(p.larvaTailS,p.larvaTailT,...
    p.larvaTailUY,p.larvaTailS(1,:),nT);
clear nT

% Set tail acceleration
% -----------------------------------------------------
nT = state.t - floor(state.t*p.kinParams.beatFreq)/p.kinParams.beatFreq;
state.tailax  = interp2(p.larvaTailS,p.larvaTailT,...
    p.larvaTailAX,p.larvaTailS(1,:),nT);
state.tailay  = interp2(p.larvaTailS,p.larvaTailT,...
    p.larvaTailAY,p.larvaTailS(1,:),nT);
clear nT

% Calculate new center of mass/volume
% -----------------------------------------------------
[Mass.net, Mass.x, Mass.y, Mass.z] = bodyMass(p.trunkParams,...
    p.ocellParams,p.finParams,p.meatParams,state.tailrx,state.tailry);

[Vol.net, Vol.x, Vol.y, Vol.z] = bodyVolume(p.trunkParams,...
    p.ocellParams,p.finParams,p.meatParams,state.tailrx,state.tailry);

I = inertiaTensor(p.trunkParams,p.ocellParams,p.finParams,...
    p.meatParams,state.tailrx,state.tailry, Mass);

state.I  = I;
state.vol = Vol.net;
state.mass = Mass.net;
state.CV = [Vol.x;Vol.y;Vol.z];
state.CM = [Mass.x;Mass.y;Mass.z];

clear Mass Vol I

% Calculate Components of Force (body frame)
% -----------------------------------------------------
[F_Tail,M_Tail]    = CalcTailForce(p,state);
[F_Trunk,M_Trunk]  = CalcTrunkForce(p,state);
[F_Buoy,M_Buoy]    = CalcBouyForce(p,state);
[F_Grav,M_Grav]    = CalcGravForce(p,state);

F_Net = F_Tail + F_Trunk + F_Buoy + F_Grav;
M_Net = M_Tail + M_Trunk + M_Buoy + M_Grav;

% Convert to inertial frame of reference
% -----------------------------------------------------
F_Net = BodyToInertial(F_Net,state.xi,state.theta,state.phi);
M_Net = BodyToInertial(M_Net,state.xi,state.theta,state.phi);

% Calculate Rotation Acceleration
% -----------------------------------------------------
dwdt = [inv(state.I)*(M_Net-cross(state.w,state.I*state.w))];

% Calculate Change in Euler Angles
% -----------------------------------------------------
dxidt     = state.w(1)+state.w(3)*cos(state.xi)*tan(state.theta)+...
    state.w(2)*sin(state.xi)*tan(state.theta);
dthetadt  = state.w(2)*cos(state.xi)-state.w(3)*sin(state.xi);
dphidt    = sec(state.theta)*(state.w(3)*cos(state.xi)+state.w(2)*sin(state.xi));

% "f" gives the derivative of the respective state variable
% -----------------------------------------------------
f      =  zeros(12,1);

f(1)   =  state.u(1);              % first derivative of position in x
f(2)   =  state.u(2);              % first derivative of position in y
f(3)   =  state.u(3);              % first derivative of position in z

f(4)   =  F_Net(1)/state.mass;     % acceleration in x
f(5)   =  F_Net(2)/state.mass;     % acceleration in y
f(6)   =  F_Net(3)/state.mass;     % acceleration in z

f(7)   =  dwdt(1);                 % rotational accel. in p
f(8)   =  dwdt(2);                 % rotational accel. in q
f(9)   =  dwdt(3);                 % rotational accel. in r

f(10)  =  dxidt;                   % changes in the Euler angles
f(11)  =  dthetadt;
f(12)  =  dphidt;


% Display function
% -----------------------------------------------------
function status = DisplayData(t,y,flag,p)
if strcmp(flag,'done')
    status = odeplot([],[],flag);
else
    status = odeplot(t,y(1:3,:),flag);
end

% Skip initializing call
if (isempty(t) | t(1) < eps) return; end;

% Set light intensity
% -----------------------------------------------------
global gBright gBrightT
for i=1:length(t)
    light(1,1) = interp1(p.lightT,p.light(1,:),t(i));
    light(2,1) = interp1(p.lightT,p.light(2,:),t(i));
    light(3,1) = interp1(p.lightT,p.light(3,:),t(i));
    light = InertialToBody(light,y(10),y(11),y(12));
    bright = -dot(light,p.ocellParams.direction);
    
    gBright = [gBright bright];
    gBrightT = [gBrightT t(i)];
end

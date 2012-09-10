function [F,M] = CalcTrunkForce(p,state)
  Re = 30;
  Cd = 24/Re + 6/(1+Re.^.5) + 0.4; %from Vogel, p.334
  Cd=20;
 
  dv = 0.5*(max(p.trunkParams.dorsal) + max(p.trunkParams.ventral));
  rl = 0.5*(max(p.trunkParams.right) + max(p.trunkParams.left));
  ap = 0.5*(min(p.trunkParams.s) + max(p.trunkParams.s));
 
  % is this projected or what
  SAx = pi*dv*rl;
  SAy = pi*dv*ap;
  SAz = pi*rl*ap;

  WBody = InertialToBody(state.w,state.xi,state.theta,state.phi);
  UBody = InertialToBody(state.u,state.xi,state.theta,state.phi);
  UBody = UBody + cross(WBody,p.larvaTrunkCV-state.CM);
 
  F(1,1)  = -0.5*p.fluidDensity*SAx*Cd*abs(UBody(1)).*UBody(1);
  F(2,1)  = -0.5*p.fluidDensity*SAy*Cd*abs(UBody(2)).*UBody(2);
  F(3,1)  = -0.5*p.fluidDensity*SAz*Cd*abs(UBody(3)).*UBody(3);

  M       = cross(p.larvaTrunkCV-state.CM,F);

function [F,M] = CalcBouyForce(p,state)
  g = InertialToBody(p.gravity,state.xi,state.theta,state.phi);
  F = -p.fluidDensity*state.vol*p.gravity;
  r = state.CV - state.CM; 
  M = cross(r,F);
  

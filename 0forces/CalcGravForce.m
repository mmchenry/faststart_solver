function [F,M] = CalcGravForce(p,state)
  g = InertialToBody(p.gravity,state.xi,state.theta,state.phi);
  F = g*state.mass;
  M = 0;

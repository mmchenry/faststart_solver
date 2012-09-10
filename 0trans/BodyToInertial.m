function out = BodyToInertial(in,xi,theta,phi);

  TM(1,1) = cos(theta)*cos(phi);
  TM(1,2) = cos(theta)*sin(phi);
  TM(1,3) = -sin(theta);
  
  TM(2,1) = sin(xi)*sin(theta)*cos(phi)-cos(xi)*sin(phi);
  TM(2,2) = sin(xi)*sin(theta)*sin(phi)+cos(xi)*cos(phi);
  TM(2,3) = sin(xi)*cos(theta);
  
  TM(3,1) = cos(xi)*sin(theta)*cos(phi)+sin(xi)*sin(phi);
  TM(3,2) = cos(xi)*sin(theta)*sin(phi)-sin(xi)*cos(phi);
  TM(3,3) = cos(xi)*cos(theta);  

  out = TM'*in;

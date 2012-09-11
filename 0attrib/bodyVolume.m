function [V,x,y,z] = bodyVolume(ocellParams,finParams,meatParams,tailX,tailY)
% This function calculates the center of volume,
% for the given morphological parameters
%
% Input data structure
%
%  ocellParams
%    -> antPost - scalar
%    -> leftRight - scalar
%    -> dorsoVent - scalar
%    -> radius - scalar
%    -> density - scalar
%
%  finParams
%    -> s - row vector (same length as tailX & tailY)
%    -> height - row vector
%    -> depth - row vector
%    -> width - row vector
%    -> density - scalar
%
%  meatParams
%    -> s - row vector (same length as tailX & tailY)
%    -> radius - row vector
%    -> density -scalar
%
%  tailX, tailY - row vector (same length as meatParams.s &
%     finParams.s), all start from x = -inf

%
% Output data structure
%  M - mass of the body
%  x,y,z - position of the center of mass

% TODO: Switch circular ocellus with elliptical swim bladder


%% Calculation of total mass

dorsalFinVol = trapz(finParams.s,...
    finParams.width.*(finParams.height-meatParams.radius));

ventralFinVol = trapz(finParams.s,...
    finParams.width.*(finParams.depth-meatParams.radius));

tailMeatVol = trapz(meatParams.s,...
    pi*meatParams.radius.^2);

ocellusVol = (4*pi/3)*ocellParams.radius^3;

V = dorsalFinVol + ventralFinVol + tailMeatVol + ocellusVol;


%% Calculation of the x position of the center of mass

dorsalFinX = trapz(finParams.s,...
    finParams.width.*(finParams.height-meatParams.radius).*tailX);

ventralFinX = trapz(finParams.s,...
    finParams.width.*(finParams.depth-meatParams.radius).*tailX);

tailMeatX = trapz(meatParams.s,...
    pi*meatParams.radius.^2.*tailX);

ocellusX = (4*pi/3)*ocellParams.radius^3*ocellParams.antPost;

x = (dorsalFinX + ventralFinX + tailMeatX + ocellusX)/V;

%% Calculation of the y position of the center of mass

dorsalFinY = trapz(finParams.s,...
    [finParams.width.*(finParams.height-meatParams.radius).*tailY]);

ventralFinY= trapz(finParams.s,...
    [finParams.width.*(finParams.depth-meatParams.radius).*tailY]);

tailMeatY = trapz(meatParams.s,...
    pi*meatParams.radius.^2.*tailY);

% By assumption of frontal symmetry

ocellusY = (4*pi/3)*ocellParams.radius^3*ocellParams.leftRight;

y = (dorsalFinY + ventralFinY + tailMeatY + ocellusY)/V;


%% Calculation of the z position of the center of mass

dorsalFinZ = trapz(finParams.s,...
    finParams.width.*(finParams.height-meatParams.radius).*...
    ((finParams.height-meatParams.radius)/2+meatParams.radius));

ventralFinZ = -trapz(finParams.s,...
    finParams.width.*(finParams.depth-meatParams.radius).*...
    0.5.*(finParams.depth+meatParams.radius));

tailMeatZ = 0;      % By assumption of circular cross-section

ocellusZ = (4*pi/3)*ocellParams.radius^3*ocellParams.dorsoVent;

z = (dorsalFinZ + ventralFinZ + tailMeatZ + ocellusZ)/V;

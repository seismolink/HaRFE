function [r]=make_rotator(strike,dip,nlay)
% Build list of rotator matrices from a list of interface strikes
% and dips. See make_rotator.m for more info.

% Topmost rotator is the identity matrix

r = zeros(3,3,nlay);
r(:,:,1) = eye(3);


% Build other rotators
layer = 2:nlay;

r(1,1,layer)=cos(strike(layer));
r(2,1,layer)=sin(strike(layer));
r(3,1,layer)=0;

r(1,2,layer)=-cos(dip(layer)).*sin(strike(layer));
r(2,2,layer)=cos(dip(layer)).*cos(strike(layer));
r(3,2,layer)=sin(dip(layer));

r(1,3,layer)=sin(dip(layer)).*sin(strike(layer));
r(2,3,layer)=-sin(dip(layer)).*cos(strike(layer));
r(3,3,layer)=cos(dip(layer));

end
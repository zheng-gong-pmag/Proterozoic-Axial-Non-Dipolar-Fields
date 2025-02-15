function [Dec,Inc,R]=Cart2Dir(x,y,z)
% Converts a direction in Cartesian coordinates into declinations and inclination.

rad=pi/180;

Dec=mod(atan2(y,x)./rad,360);
Inc=atan2(z,sqrt(x.*x+y.*y))./rad;
R=sqrt(x.*x+y.*y+z.*z);

end
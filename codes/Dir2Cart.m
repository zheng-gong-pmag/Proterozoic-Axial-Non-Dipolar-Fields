function [x,y,z]=Dir2Cart(Dec,Inc)
% Converts a list or array of vector directions in degrees (declination, inclination)
% to an array of the direction in Cartesian coordinates (x,y,z).
rad=pi/180;

x=cos(rad.*Inc).*cos(rad.*Dec);
y=cos(rad.*Inc).*sin(rad.*Dec);
z=sin(rad.*Inc);

end
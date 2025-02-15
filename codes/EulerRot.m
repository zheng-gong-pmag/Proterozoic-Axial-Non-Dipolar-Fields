% Euler rotation
% by Zheng Gong, Dec 7, 2019

function [Rot_Plat,Rot_Plon] = EulerRot(Plat,Plon,EulerLat,EulerLon,EulerAng)

rad=pi/180;

% Convert Lat/Lon to Cartesian
x=cos(rad.*Plat).*cos(rad.*Plon);
y=cos(rad.*Plat).*sin(rad.*Plon);
z=sin(rad.*Plat);

Euler_x=cos(rad.*EulerLat).*cos(rad.*EulerLon);
Euler_y=cos(rad.*EulerLat).*sin(rad.*EulerLon);
Euler_z=sin(rad.*EulerLat);

% Compute rotation matrix elements
Rot11=Euler_x.*Euler_x.*(1-cos(rad.*EulerAng))+cos(rad.*EulerAng);
Rot12=Euler_x.*Euler_y.*(1-cos(rad.*EulerAng))-Euler_z.*sin(rad.*EulerAng);
Rot13=Euler_x.*Euler_z.*(1-cos(rad.*EulerAng))+Euler_y.*sin(rad.*EulerAng);

Rot21=Euler_y.*Euler_x.*(1-cos(rad.*EulerAng))+Euler_z.*sin(rad.*EulerAng);
Rot22=Euler_y.*Euler_y.*(1-cos(rad.*EulerAng))+cos(rad.*EulerAng);
Rot23=Euler_y.*Euler_z.*(1-cos(rad.*EulerAng))-Euler_x.*sin(rad.*EulerAng);

Rot31=Euler_z.*Euler_x.*(1-cos(rad.*EulerAng))-Euler_y.*sin(rad.*EulerAng);
Rot32=Euler_z.*Euler_y.*(1-cos(rad.*EulerAng))+Euler_x.*sin(rad.*EulerAng);
Rot33=Euler_z.*Euler_z.*(1-cos(rad.*EulerAng))+cos(rad.*EulerAng);

% Apply rotation
Rot_x=Rot11.*x+Rot12.*y+Rot13.*z;
Rot_y=Rot21.*x+Rot22.*y+Rot23.*z;
Rot_z=Rot31.*x+Rot32.*y+Rot33.*z;

% Convert Cartesian to Lat/Lon
Rot_Plat=asin(Rot_z)./rad;
Rot_Plon=atan(Rot_y./Rot_x)./rad;

Rot_Plon(Rot_x >= 0) = mod(Rot_Plon(Rot_x >= 0),360);
Rot_Plon(Rot_x < 0) = mod(Rot_Plon(Rot_x < 0)+180,360);

end
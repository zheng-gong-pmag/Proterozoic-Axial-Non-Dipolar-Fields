function [theta]=AngDiff(Dec1,Inc1,Dec2,Inc2)

rad=pi/180;

xx=cos(Inc1.*rad).*cos(Dec1.*rad).*cos(Inc2.*rad).*cos(Dec2.*rad);
yy=cos(Inc1.*rad).*sin(Dec1.*rad).*cos(Inc2.*rad).*sin(Dec2.*rad);
zz=sin(Inc1.*rad).*sin(Inc2.*rad);

theta=acos(xx+yy+zz)./rad;
theta=real(theta);

end
function [Inclination]=getInclinationG2G3(Latitude,G2,G3)


phi=90-Latitude;
% I_num=2.*cosd(phi)+1.5.*G2.*(3*cosd(phi).^2-1)+2.*G3.*(5.*cosd(phi).^3-3.*cosd(phi));
% I_den=sind(phi)+G2.*(3.*sind(phi).*cosd(phi))+1.5.*G3.*(5.*sind(phi).*cosd(phi).^2-sind(phi));
% Inclination=atand(I_num./abs(I_den));


% pre-calcuate the params that will be used repeatedly to reduce computing time

cos_phi=cosd(phi);
sin_phi=sind(phi);

cos_phi_2=cos_phi.^2;
sin_phi_cos_phi=sin_phi.*cos_phi;
sin_phi_cos_phi_2=sin_phi.*cos_phi_2;

I_num=2*cos_phi+1.5.*G2*(3*cos_phi_2-1)+2.*G3*(5*cos_phi_2.*cos_phi-3*cos_phi);
I_den=sin_phi+G2.*(3*sin_phi_cos_phi)+1.5.*G3*(5*sin_phi_cos_phi_2-sin_phi);

Inclination=atand(I_num./abs(I_den));
% Inclination=atand(I_num./I_den);

end
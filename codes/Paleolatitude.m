function [Paleolatitude,PLatGAD]=Paleolatitude(Inclination, G_2, G_3)

rad=pi/180;

PLatGAD=atan(0.5.*tan(rad.*Inclination))./rad;

for phi_deg=0:0.01:180
    
    phi=phi_deg.*rad;
    I_num=2.*cos(phi)+1.5.*G_2.*(3*cos(phi).^2-1)+2.*G_3.*(5.*cos(phi).^3-3.*cos(phi));
    I_den=sin(phi)+G_2.*(3.*sin(phi).*cos(phi))+1.5.*G_3.*(5.*sin(phi).*cos(phi).^2-sin(phi));
    Calculated_Inclination=atan(I_num./I_den)./rad;
    Paleolatitude=90-(phi./rad);
    if abs(Inclination-Calculated_Inclination)<=0.1
        break
    end
    
end

end
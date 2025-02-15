function [Plat,Plon,A95,dp,dm,paleolat]=Dir2VGP(Slat,Slon,Dm,Im,a95)

% VGP calculation from site-level direction

rad=pi/180;           % degree to radian

Slat=Slat*rad;        % site latitude  N
Slon=Slon*rad;        % site longitude E
Dm=Dm*rad;            % site-mean declination
Im=Im*rad;            % site-mean inclination

% calculations
p =	atan(tan(Im)*0.5);
Plat = 90-acos(sin(Slat).*sin(p)+cos(Slat).*cos(Dm).*cos(p))/rad;
beta = asin(cos(p).*sin(Dm)./cos(Plat*rad));

A = sin(p);
B = sin(Slat).*sin(Plat.*rad);

Plon=nan(length(Dm),1);

for i=1:length(Dm)

    if	A(i) < B(i)
    	Plon(i,1) = mod(((Slon(i)+pi-beta(i))/rad),360);
    else
    	Plon(i,1) = mod(((Slon(i)+beta(i))/rad),360);
    end
    
end

dp = a95.*((1+3*sin(p).*sin(p))/2);
dm = a95.* cos(p)./cos(Im);
A95 = sqrt(dp.*dm);
paleolat = p/rad;                     %	paleo-latitude

end
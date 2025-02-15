function [PoleLat, PoleLong, PaleoSiteIncl] = Pole(Paleolat, SiteLat, SiteLong, SiteDecl)

rad=pi/180;
N=length(Paleolat);

PaleoSiteIncl=atan(2*tan(Paleolat.*rad))./rad;

% Calculate Pole latitude based on great circle trigonometry, uses North Pole as apex
PoleLat=90-(acos(cos(rad.*(90-SiteLat)).*cos(rad.*(90-Paleolat))+sin(rad.*(90-SiteLat)).*sin(rad.*(90-Paleolat)).*cos(rad.*SiteDecl)))./rad;

% Butler's beta parameter (angle made between great circle paths).
beta=asin(sin(rad.*SiteDecl).*sin(rad.*(90-Paleolat))./sin(rad.*(90-PoleLat)))./rad;

% Ambiguity test for longitude
ambiguity=cos(rad.*(90-Paleolat))-sin(rad.*PoleLat).*sin(rad.*SiteLat);

PoleLong=nan(N,1);

% Calculate Pole longitude in and angle between 0 and 360 from East
for i=1:N
    
    if ambiguity(i,1)<0
        PoleLong(i,1)=mod(SiteLong(i,1)+180-beta(i,1),360);
    else
        PoleLong(i,1)=mod(SiteLong(i,1)+beta(i,1),360);
    end

    if PoleLat(i,1)==90
        PoleLong(i,1)=0;
    end

end

% % Put the poles in a single polarity for comparison
% for i=1:N
%     if PoleLong(i,1)>270
%         PoleLat(i,1)=-1*PoleLat(i,1);
%         PoleLong(i,1)=mod(PoleLong(i,1)+180,360);
%     end
%     if PoleLong(i,1)<90
%         PoleLat(i,1)=-1*PoleLat(i,1);
%         PoleLong(i,1)=mod(PoleLong(i,1)+180,360);
%     end
% end

end
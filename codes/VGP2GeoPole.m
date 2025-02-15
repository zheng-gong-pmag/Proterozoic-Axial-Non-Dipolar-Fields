% rotate virtual geomagnetic pole to north or south pole

function [EulerLat,EulerLon,EulerAng]=VGP2GeoPole(PoleLat,PoleLon,PoleType)

EulerLat=0.*PoleLat;

condition1=strcmpi(PoleType,{'North','N'});
condition2=strcmpi(PoleType,{'South','S'});

if sum(condition1)==1

    EulerLon=mod(PoleLon-90,360);
    EulerAng=90-PoleLat;

elseif sum(condition2)==1

    EulerLon=mod(PoleLon+90,360);
    EulerAng=90+PoleLat;

end

end
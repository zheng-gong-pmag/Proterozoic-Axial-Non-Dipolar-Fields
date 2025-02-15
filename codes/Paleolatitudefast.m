function [Paleolatitude]=Paleolatitudefast(Inclination,Latitude,G2,G3)

latitudesteps=(-90:0.1:90);

modeled_Inclination=getInclinationG2G3(latitudesteps,G2,G3);

% find 10 most likely intersections
[~,candidate_id]=mink(abs(Inclination-modeled_Inclination),10);

Paleolatitude_candidates=latitudesteps(candidate_id);

% find the most likely intersection (to solve the non-unique solution issue)
[~,id]=min(abs(Paleolatitude_candidates-Latitude));

Paleolatitude=Paleolatitude_candidates(id);

end
function [A95]=KN2A95(K,N)

rad=pi/180;

R=N-((N-1)./K);
A95=acos(1-((N-R)./R).*(20^(1/(N-1))-1))./rad;

end
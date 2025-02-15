function [c]=Val2Clr(aid,amax,amin,clrmap)

cid=round(255.*(aid-amin)./(amax-amin))+1;
c=clrmap(cid,:);

end
function [p]=zeroinCI2(data)

n=length(data);

data_sort=sort(data);

medianvalue=median(data_sort);

p=1;
h=zeros(100,1);

for i=1:1:100

    CIpercentage=1-i/100;

    minusvalue=data_sort(ceil(n*(0.5-CIpercentage/2)));

    plusvalue=data_sort(floor(n*(0.5+CIpercentage/2)));

    if minusvalue*plusvalue>0

        h(i,1)=1;
        p=1-CIpercentage;

        break
    end

end

end
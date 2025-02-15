function [h,medianvalue,plusvalue,minusvalue]=zeroinCI(data,kwargs)

arguments
    data;
    kwargs.sigma = 2;
end


if kwargs.sigma == 2

    CIpercentage=0.95;

elseif kwargs.sigma == 1

    CIpercentage=0.68;

elseif kwargs.sigma == 3

    CIpercentage=0.997;

end


n=length(data);

data_sort=sort(data);

medianvalue=median(data_sort);

minusvalue=data_sort(ceil(n*(0.5-CIpercentage/2)));

plusvalue=data_sort(floor(n*(0.5+CIpercentage/2)));

h=1;

if minusvalue*plusvalue<=0

    h=0;

end

end
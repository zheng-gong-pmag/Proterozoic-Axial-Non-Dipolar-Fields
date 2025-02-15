function [Dec_f,Inc_f,Group]=FlipDir(Dec,Inc,kwargs)

arguments
    Dec;
    Inc;
    kwargs.positive {MustBeBoolean(kwargs.positive)} = true;
end

N=length(Dec);

anglediff=nan(N,1);
Dec_f=nan(N,1);
Inc_f=nan(N,1);
Group=nan(N,1);

if kwargs.positive
    id=(Inc>=0);
else
    id=(Inc<0);
end

for i=1:N

    anglediff(i,1)=AngDiff(median(Dec(id)),median(Inc(id)),Dec(i),Inc(i));

    if anglediff(i,1)>90
        Dec_f(i,1)=mod(Dec(i)+180,360);
        Inc_f(i,1)=-Inc(i);
        Group(i,1)=0;
    else
        Dec_f(i,1)=Dec(i);
        Inc_f(i,1)=Inc(i);
        Group(i,1)=1;
    end

end

end

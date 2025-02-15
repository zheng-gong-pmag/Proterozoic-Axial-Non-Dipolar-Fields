function [tau]=getElongation(Dec,Inc)

[x1,x2,x3]=Dir2Cart(Dec,Inc);

T11=sum(x1.*x1);  T12=sum(x1.*x2);  T13=sum(x1.*x3);
T21=sum(x2.*x1);  T22=sum(x2.*x2);  T23=sum(x2.*x3);
T31=sum(x3.*x1);  T32=sum(x3.*x2);  T33=sum(x3.*x3);

T=[T11,T12,T13;T21,T22,T23;T31,T32,T33];

[~,D]=eig(T);

[Ddiag]=sort(diag(D),'descend');

e2=Ddiag(2);
e3=Ddiag(3);
tau=e2./e3;

end

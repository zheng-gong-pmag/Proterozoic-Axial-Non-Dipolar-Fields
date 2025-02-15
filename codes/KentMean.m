function [DecM,IncM,A95M,KM,ZDec,ZInc,EDec,EInc,Zeta,Eta,Kappa,Beta] = KentMean(Dec,Inc,NN,distribution_95)
% Gets Kent parameters for data.

% Parameters:
%   data : nested pairs of [Dec, Inc]
%   NN : normalization
%       Number of data for Kent ellipse
%       NN is 1 for Kent ellipses of bootstrapped mean directions
%   distribution_95 : the default behavior (distribution_95=false) is for
%       the function to return the confidence region for the mean direction.
%       If distribution_95=true, what will be returned are the parameters
%       associated with the region containing 95% of the directions.

% Returns:
%   kpars : dictionary keys
%       DecM  : mean declination
%       IncM  : mean inclination
%       A95M  : Fisher A95
%       KM    : Fisher precision parameter
%       Zeta  : major ellipse
%       ZDec  : declination of major ellipse axis
%       ZInc  : inclination of major ellipse axis
%       Eta   : minor ellipse
%       EDec  : declination of minor ellipse axis
%       EInc  : inclination of minor ellipse axis
%       Kappa : Kent concentration parameter
%       Beta  : Kent ovalness parameter

N = length(Dec);

if N < 2
    return;
end

% get fisher mean and convert to co-inclination (theta)/dec (phi) in radians
[DecM,IncM,A95M,KM,R] = FisherMean(Dec,Inc);
pbar = DecM * pi / 180;
tbar = (90 - IncM) * pi / 180;

% initialize matrices
H = [cos(tbar) * cos(pbar), -sin(pbar), sin(tbar) * cos(pbar);
    cos(tbar) * sin(pbar), cos(pbar), sin(pbar) * sin(tbar);
    -sin(tbar), 0, cos(tbar)];

% get cartesian coordinates of data
X = nan(N,3);
gam=zeros(3);

for i = 1:N
    [X(i,1),X(i,2),X(i,3)] = Dir2Cart(Dec(i),Inc(i));
end

% put in T matrix
T = Tmatrix(X);
for i = 1:3
    for j = 1:3
        T(i,j) = T(i,j) / NN;
    end
end

% compute B=H'TH
w = zeros(3);
b = zeros(3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            w(i,j) = w(i,j) + T(i,k) * H(k,j);
        end
    end
end

for i = 1:3
    for j = 1:3
        for k = 1:3
            b(i,j) = b(i,j) + H(k,i) * w(k,j);
        end
    end
end

% get the upper 2x2 matrix of B as BL
BL = b([1, 2], [1, 2]);

% compute the eigenvalues of BL
BL_ei = eig(BL);

% choose a rotation w about North pole to diagonalize upper part of B
psi = 0.5 * atan(2 * b(1,2) / (b(1,1) - b(2,2)));
w = [cos(psi), -sin(psi), 0;
    sin(psi), cos(psi), 0;
    0, 0, 1];

for i=1:3
    for j=1:3
        gamtmp=0;
        for k=1:3
            gamtmp=gamtmp+H(i,k).*w(k,j);
            
        end
        gam(i,j)=gamtmp;
    end
end
      
xg = zeros(N, 3);

for i = 1:N
    xg(i, :) = [0, 0, 0];
    for k = 1:3
        xgtmp = 0;
        for j = 1:3
            xgtmp = xgtmp + gam(j, k) * X(i, j);
        end
        xg(i, k) = xgtmp;
    end
end

% compute asymptotic ellipse parameters
xmu = 0;
sigma1 = 0;
sigma2 = 0;

for i = 1:N
    xmu = xmu + xg(i,3);
    sigma1 = sigma1 + xg(i,1)^2;
    sigma2 = sigma2 + xg(i,2)^2;
end

xmu = xmu / N;
sigma1 = sigma1 / N;
sigma2 = sigma2 / N;

if sum(strcmpi(distribution_95,{'false','f'})) == 1
    g = -2.0 * log(0.05) / (NN * xmu^2);
elseif sum(strcmpi(distribution_95,{'true','t'})) == 1
    g = -2.0 * log(0.05) / (xmu^2);
end


if sqrt(sigma1 * g) < 1
    eta = asin(sqrt(sigma1 * g));
end

if sqrt(sigma2 * g) < 1
    zeta = asin(sqrt(sigma2 * g));
end

if sqrt(sigma1 * g) >= 1
    eta = pi / 2.0;
end

if sqrt(sigma2 * g) >= 1
    zeta = pi / 2.0;
end

if eta > zeta
    
    eta_output=zeta;
    zeta_output=eta;

    [EDec,EInc,~] = Cart2Dir(gam(1,2), gam(2,2), gam(3,2));
    [ZDec,ZInc,~] = Cart2Dir(gam(1,1), gam(2,1), gam(3,1));

else
    
    eta_output=eta;
    zeta_output=zeta;

    [ZDec,ZInc,~] = Cart2Dir(gam(1,2), gam(2,2), gam(3,2));
    [EDec,EInc,~] = Cart2Dir(gam(1,1), gam(2,1), gam(3,1));

end

if ZInc < 0
    ZInc = -ZInc;
    ZDec = mod(ZDec + 180, 360);
end
if EInc < 0
    EInc = -EInc;
    EDec = mod(EDec + 180, 360);
end
Zeta = zeta_output * 180 / pi;
Eta = eta_output * 180 / pi;
R1 = R / N;
R2 = abs(BL_ei(1) - BL_ei(2));

Kappa=1/(2-2*R1-R2)+1/(2-2*R1+R2);
Beta=0.5*(1/(2-2*R1-R2)-1/(2-2*R1+R2));

end

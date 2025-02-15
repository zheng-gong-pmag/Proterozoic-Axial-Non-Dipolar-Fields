function result = kent_mle(x)
    % x is the data in Euclidean coordinates
    tic;
    n = size(x, 1);  % sample size
    xbar = mean(x);  % mean vector
    S = x' * x / n;
    xbar = xbar / sqrt(sum(xbar.^2)); % mean direction
    u = [acos(xbar(1)), atan2(xbar(3), xbar(2)) + pi * (xbar(2) < 0)] / (2 * pi);
    % u is the mean vector to latitude and longitude
    theta = u(1);
    phi = u(2);
    costheta = cos(theta);
    sintheta = sin(theta);
    cosphi = cos(phi);
    sinphi = sin(phi);
    H = [costheta, sintheta * cosphi, sintheta * sinphi; ...
        -sintheta, costheta * cosphi, costheta * sinphi; ...
        0, -sinphi, cosphi];
    B = H' * S * H;
    psi = 0.5 * atan(2 * B(2, 3) / (B(2, 2) - B(3, 3)));
    K = [1, 0, 0; 0, cos(psi), sin(psi); 0, -sin(psi), cos(psi)];
    G = H * K;  % The G matrix Kent describes, the A in our notation
    lam = eig(B(2:end, 2:end));
    
    % the next function will be used to estimate the kappa and beta
    xg = x * G;
    xg1 = sum(xg(:, 1));
    a = sum(xg(:, 2:3).^2);
    xg2 = a(1);
    xg3 = a(2);
    
    % Maximization w.r.t. to k and b
    mle = @(para) mle_func(para, n, xg1, xg2, xg3);
%     ini = vmf_mle(x);
%     ini = [ini, ini / 2.1];  % initial values for kappa and beta
    [Dec,Inc,~]=Cart2Dir(x(:,1),x(:,2),x(:,3));
    [~,~,~,~,~,~,~,~,~,~,Kappa,~] = KentMean(Dec,Inc,length(Dec), 'f');
    ini=[Kappa,Kappa/2.1];


    qa = fminunc(mle, ini);
    para = qa;
    k = para(1);
    b = para(2);  % the estimated parameters
    gam = [0, k, 0];
    lam = [0, -b, b];
    ckb = fb_saddle(gam, lam);
    % the line below calculates the log-likelihood
    l = -n * ckb + k * xg1 + b * (xg2 - xg3);
    param = [k, b, psi];
    runtime = toc;
    result.G = G;
    result.param = param;
    result.logcon = ckb;
    result.loglik = l;
    result.runtime = runtime;
end

function val = mle_func(para, n, xg1, xg2, xg3)
    k = para(1);
    b = para(2);
    gam = [0, k, 0];
    lam = [0, -b, b];
    ckb = fb_saddle(gam, lam);
    ckb=ckb(3);
    val = n * ckb - k * xg1 - b * (xg2 - xg3);
end

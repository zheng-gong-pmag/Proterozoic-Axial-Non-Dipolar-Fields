function logcon = fb_saddle(gam, lam)
    % Saddlepoint approximations of the Fisher-Bingham distributions
    % Tsagris Michail 02/2014
    % References: Kume Alfred and Wood Andrew T.A. (2005)
    % Saddlepoint approximations for the Bingham and Fisher-Bingham normalizing constants (Biometrika)

    lam = sort(lam);  % sorts the eigenvalues of the Bingham part
    mina = min(lam);
    if mina <= 0
        aaa = abs(mina) + 1;
        lam = lam + aaa;  % makes all the lambdas positive and greater than zero
    end
    p = length(gam);  % dimensionality of the distribution
    para = [gam, lam];  % the parameters of the Fisher-Bingham

    % Saddlepoint equation
    saddle_equat = @(ta, para) sum(0.5 ./ (lam - ta) + 0.25 * (gam.^2 ./ (lam - ta).^2)) - 1;

    low = lam(1) - 0.25 * p - 0.5 * sqrt(0.25 * p^2 + p * max(gam)^2);  % lower bound
    up = lam(1) - 0.25 - 0.5 * sqrt(0.25 + min(gam)^2);  % not the exact upper bound but a bit higher
    tau = fzero(@(ta) saddle_equat(ta, para), [low, up]);  % tau which solves the saddlepoint equation

    % Derivatives of the cumulant generating function
    kfb = @(j, gam, lam, ta) sum(0.5 * factorial(j-1) ./ (lam - ta).^j + 0.25 * factorial(j) * gam.^2 ./ (lam - ta).^(j + 1));

    rho3 = kfb(3, gam, lam, tau) / kfb(2, gam, lam, tau)^1.5;
    rho4 = kfb(4, gam, lam, tau) / kfb(2, gam, lam, tau)^2;
    Ta = rho4 / 8 - 5 / 24 * rho3^2;
    c1 = 0.5 * log(2) + 0.5 * (p - 1) * log(pi) - 0.5 * log(kfb(2, gam, lam, tau)) - 0.5 * sum(log(lam - tau)) - tau + 0.25 * sum(gam.^2 ./ (lam - tau));
    c2 = c1 + log1p(Ta);
    c3 = c1 + Ta;

    % If there were negative lambdas, adjust the log constants
    if mina <= 0
        c1 = c1 + aaa;
        c2 = c2 + aaa;
        c3 = c3 + aaa;
    end

    logcon = [c1, c2, c3];
    names = {'first order', 'second order', 'third order'};
    logcon = array2table(logcon, 'VariableNames', names);
end

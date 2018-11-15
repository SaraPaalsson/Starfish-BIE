function self_contained_est_demo()
% Demo of error estimates for Laplace double layer potential
%
% (C) 2018 Ludvig af Klinteberg
    
    % resolution of plot grid
    gridpoints = 150;
    
    % Parametrization of panel
    % k=Curving of panel
    k = 0.4; 
    g = @(t) t + 1i*k*t.^2;
    gp = @(t) 1 + 1i*2*k*t;    
    n = @(t) 1i*gp(t)./abs(gp(t));
    
    % Density to use
    density_g = @(t) real(g(t)).*imag(g(t)); % sigma(t) = x(t)*y(t)
    density_t = @(t) t; % sigma(t) = t
    density_1 = @(t) t.^0; % sigma(t) = 1
    
    
    density = density_g; % try g/t/1
    
    % Laplace double layer kernel
    kernel = @(x, y, ny) real( ny ./ (x-y) );
    integrand = @(t,z) kernel(z, g(t), n(t)) ...
        .* density(t) .* abs(gp(t));
    
    
    % Gauss-Legendre quadrature nodes and weights
    [t16,w16] = lgwt(16,-1,1);
    d16 = density(t16);
    
    % Create Legendre expansion of panel (in rescaled coordinates)
    za = g(-1);
    zb = g(1);
    L = legendre_matrix(16);
    coeffs = L*rotate_and_scale(za, zb, g(t16));
    % Create expansion of density on panel
    den_coeffs = L*d16;
    
    % Create a grid for evaluation
    x = linspace(-2.0, 2.0, gridpoints);
    y = linspace(-1.5, 1, gridpoints);
    [X,Y] = meshgrid(x,y);
    Z = X + 1i*Y;
    
    % Compute values on grid
    [I, Q16, est] = deal(zeros(size(Z)));
    n = 16;
    for i=1:numel(Z)
        z = Z(i);
        % === Reference
        % Kill warnings about resolution not being enough
        warning off
        I(i) = integral(@(t) integrand(t, z), -1, 1, ...
                        'abstol', eps, 'reltol', eps);
        warning on
        % === Direct
        Q16(i) = sum(integrand(t16, z) .* w16);
        % * Find root
        zr = rotate_and_scale(za, zb, z);
        tinit = zr;
        t = newton_legendre(coeffs, zr, tinit, 20, 1e-13);
        P = legendre_deriv_scalar(15, t);
        % compute density at root
        den_t = sum(P(:).*den_coeffs(:));
        sgn = sign(real(t));
        % compute estimate
        est(i) = imag(2*pi*den_t*(t + sgn*sqrt(t.^2-1)).^(-2*n-1));
    end
    
    % Compute rel. errors and normalize est
    Inorm = norm(I(:), inf);
    E16 = abs(Q16-I) / Inorm;
    est = est / Inorm;
    
    % Plot
    figure(1)
    clf
    pcolor(real(Z), imag(Z), log10(E16));    
    title('Direct quadrature (16 points)')
    
    num_levels = 7;
    levels = linspace(-14, 0, num_levels+1);
    colormap('default');
    cmap = colormap();
    idx = round(linspace(1, size(cmap,1), num_levels));
    cmap = cmap( idx, :);
    colormap(cmap);
    
    cbar = colorbar();
    set(cbar, 'Ticks', levels)
    shading interp
    caxis([min(levels), max(levels)])
    hold on
    contour(real(Z), imag(Z), log10(abs(est)), levels(2:end), 'k')
    plot(complex(g(t16)), '.-r')    
    axis equal
end

% ============= END MAIN SUB ==============

% ========== INCLUDED FUNCTIONS ==========

function zr = rotate_and_scale(za, zb, z)
% zr = rotate_and_scale(za, zb, z)
%
% zr = M(z)
% where
% M(za) = -1
% M(zb) = 1
    c_map = [ (za+zb)/(za-zb) 2/(zb-za) ];
    zr = c_map(1) + z*c_map(2);
end

function L = legendre_matrix(order)
% L = legendre_matrix(order)
% 
% L computes Legendre expansion coefficients for a function defined at the
% Gauss-Legendre quadrature nodes on [-1, 1]
%
%  c = L*f
%
% c_l = (2l+1)/2 \sum_n P_l(x_n) f_n w_n
    [x, w] = lgwt(order, -1, 1);
    P = legendre_vec(order-1, x); %P = [P_0(x) P_1(x) P_2(x) ...]
    l = 0:order-1;
    L = bsxfun(@times, P, w);
    L = bsxfun(@times, (2*l+1)/2, L);
    L = L.';
end

function P = legendre_vec(n, x)
% P = legendre_vec(n, x)
% => P(:, l+1) = P_l(x)
% 0 <= l <= n
% 
% Recurrence relation: [http://dlmf.nist.gov/18.9]
% (l+1)P_(l+1)(x)-(2l+1)xP_l(x)+lP_(l-1)(x)=0 	
    x = x(:);
    P = zeros(numel(x), n+1);
    P(:, 1) = ones(size(x)); % l=0
    P(:, 2) = x; % l=1
    for l=1:n-1
        % Compute l+1
        P(:, l+1+1) = ( (2*l+1)*x .* P(:, l+1) - l*P(:, l-1+1) ) / (l+1);
    end
end

function [P, D] = legendre_deriv_scalar(n, x)
    [P, D] = deal(complex(zeros(1, n+1)));
    P(1) = 1; % l=0
    D(1) = 0;
    P(2) = x; % l=1
    D(2) = 1;
    for l=1:n-1
        % Compute l+1
        P(l+1+1) = ( (2*l+1)*x * P(l+1) - l*P(l-1+1) ) / (l+1);
        D(l+1+1) = ( (2*l+1)*(P(l+1) + x * D(l+1)) - l*D(l-1+1) ) / (l+1);    
    end
end

function [t, relres, iter] = newton_legendre(c, zr, t0, maxiter, tol)
    t = t0;
    iter = 0;
    lmax = numel(c)-1;
    relres = tol;
    for iter=1:maxiter
        [P, D] = legendre_deriv_scalar(lmax, t);
        f = P*c - zr;
        fp = D*c;
        dt = -f / fp;
        t = t + dt;
        relres = abs(dt/t);
        if relres < tol
            return
        end
    end
end

function [x,w]=lgwt(N,a,b)
% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
    N=N-1;
    N1=N+1; N2=N+2;
    xu=linspace(-1,1,N1)';
    % Initial guess
    y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
    % Legendre-Gauss Vandermonde Matrix
    L=zeros(N1,N2);
    % Derivative of LGVM
    Lp=zeros(N1,N2);
    % Compute the zeros of the N+1 Legendre Polynomial
    % using the recursion relation and the Newton-Raphson method
    y0=2;
    % Iterate until new points are uniformly within epsilon of old points
    while max(abs(y-y0))>eps      
        L(:,1)=1;
        Lp(:,1)=0;        
        L(:,2)=y;
        Lp(:,2)=1;        
        for k=2:N1
            L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
        end        
        Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);           
        y0=y;
        y=y0-L(:,N2)./Lp;        
    end
    % Linear map from[-1,1] to [a,b]
    x=flipud((a*(1-y)+b*(1+y))/2);      
    % Compute the weights
    w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end

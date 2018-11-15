% Test estimating error around a starfish panel for laplace's equ. To
% understand Ludvig's self_contained_est_demo.m.
function estimateStokes()

clc

disp('Start estimateStokes')

gridpoints = 150;
disp(['Error computed at ' num2str(gridpoints) ' gridpoints'])

% THIS IS THE EXACT PARAMETRISATION THAT WE do not KNOW ABOUT
% === Parametrization of panel
% k=Curving of panel
k = 0.4;
g = @(t) t + 1i*k*t.^2;
gp = @(t) 1 + 1i*2*k*t;
n = @(t) 1i*gp(t)./abs(gp(t));

% === Density on panel
density_1 = @(t) t.^0; % sigma(t) = 1
density_t = @(t) t; % sigma(t) = t
density_g = @(t) real(g(t)).*imag(g(t)); % sigma(t) = x(t)*y(t)

density = density_1;

% === Laplace double layer kernel (integral 1) (obs no imag)
kernel1 = @(x, y, ty) imag(  ty./ (y-x) );
integrand1 = @(t,z) 1/(pi*1i) * kernel1(z, g(t), gp(t)) .* density(t);

% === Integral 2
kernel2 =  @(x,y,ty) ty./((y-x));
integrand2 = @(t,z) kernel2(z,g(t),gp(t)) .* conj(n(t)).^2 .* density(t);

% === Integral 3
kernel3 = @(x,y,ty) ty./((y-x).^2);
integrand3 = @(t,z) kernel3(z,g(t),gp(t)) .* conj(g(t)-z) .* density(t);


% HERE WE PARAMETRISE THE STRAIGHT LINE BETWEEN -1 AND 1
% === Gauss-Legendre quadrature nodes and weights
[t16,w16] = lgwt(16,-1,1);
d16 = density(t16);
z16 = g(t16);
zp16 = gp(t16);
n16 = 1i*zp16./abs(zp16);
% When we use this later, d16,z16,zp16 is all we know.

% === Create Legendre expansion of panel (in rescaled coordinates)
za = g(-1); % Starting point of panel
zb = g(1); % End point of panel
L = legendre_matrix(16);
coeffs = L*rotate_and_scale(za, zb, z16);

% === Create expansion of density on panel
den_coeffs = L*d16;
zp_coeff = L*zp16;
n_coeff = L*conj(n16).^2;

% === Create a grid for evaluation
x = linspace(-2.0, 2.0, gridpoints);
y = linspace(-1.5, 1, gridpoints);
[X,Y] = meshgrid(x,y);
Z = X + 1i*Y;

n = 16;

[I1, Q1, est1, I2, Q2, est2, I3, Q3, est3, Itot, Qtot, esttot] = deal(zeros(size(Z)));
for j=1:numel(Z)
    z = Z(j);
    
    % === Reference
    % Kill warnings about resolution not being enough
    warning off
    I1(j) = integral(@(t) integrand1(t, z), -1, 1, ...
        'abstol', eps, 'reltol', eps);
    I2(j) = -1/(2*pi) * conj( integral(@(t) integrand2(t, z), -1, 1, ...
        'abstol', eps, 'reltol', eps) );
    I3(j) = -1/(2*pi) * conj( integral(@(t) integrand3(t,z), -1, 1, ...
        'abstol', eps, 'reltol', eps) );
    warning on
    Itot(j) = I1(j) + I2(j) + I3(j);
    
    % === Direct
    Q1(j) = 1/(pi*1i) * sum(d16 .* kernel1(z,z16,zp16) .* w16);
    Q2(j) = -1/(2*pi) * conj( sum(d16 .* conj(n16).^2 .* kernel2(z,z16,zp16) .* w16) );
    Q3(j) = -1/(2*pi) * conj( sum(d16 .* conj(z16-z) .* kernel3(z,z16,zp16) .* w16) );
    Qtot(j) = Q1(j) + Q2(j) + Q3(j);
    
    % === Find root
    zr = rotate_and_scale(za, zb, z);
    tinit = zr;
    t0 = newton_legendre(coeffs, zr, tinit, 20, 1e-13);
    P = legendre_deriv_scalar(15, t0);
    
    % === Compute density at root
    den_t = sum(P(:).*den_coeffs(:));
    % === Compute zp at root
    zp_t = sum(P(:).*zp_coeff(:));
    % === Compute n at root
    n_t = sum(P(:).*n_coeff(:));
    % === Compute factor x(t0) - x(t0') - 1i*(y(t0) - y(t0')) = g(t0')' - g(t0)'
    gt0 = P*coeffs;
    Pconj = legendre_deriv_scalar(15,conj(t0));
    gt0conj = Pconj*coeffs;
    factor = conj(gt0conj) - conj(gt0);
    
    
    % === Compute estimate
    % qth derivative of k_n(z):
    knq = @(n, q, z) (-(2*n+1)./sqrt(z.^2-1)).^q .* 2*pi./(z+sign(real(z))*sqrt(z.^2-1)).^(2*n+1);

    est1(j) = abs( 1/(pi*1i) * imag((den_t * knq(n,0,t0))) );
    est2(j) = 1/(2*pi) * abs( conj(den_t * n_t * knq(n,0,t0)) );
    est3(j) = 1/(2*pi) * abs( conj(den_t .* knq(n,1,t0)/zp_t * factor) );
    esttot(j) = abs(est1(j) + est2(j) + est3(j));
    est1(j) = abs(est1(j)); est2(j) = abs(est2(j)); est3(j) = abs(est3(j));
        
end

%  === Compute rel. errors and normalize est
Inorm1 = norm(I1(:), inf);
E16_1 = abs(Q1-I1) / Inorm1;
est1 = est1 / Inorm1 + eps;

Inorm2 = norm(I2(:), inf);
E16_2 = abs(Q2-I2) / Inorm2;
est2 = est2 / Inorm2 + eps;

Inorm3 = norm(I3(:), inf);
E16_3 = abs(Q3-I3) / Inorm3;
est3 = est3 / Inorm3 + eps;

Inormtot = norm(Itot(:), inf);
E16_tot = abs(Qtot-Itot) / Inormtot;
esttot = esttot / Inormtot + eps;


% === Plot
figure(1)
clf
pcolor(real(Z), imag(Z), log10(E16_1));
title('Direct quadrature (16 points), integral 1')

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
contour(real(Z), imag(Z), log10(abs(est1)), levels(2:end), 'k')
plot(complex(g(t16)), '.-r')
axis equal
%-------------------------------------------------------------------------
figure(2)
clf
pcolor(real(Z), imag(Z), log10(E16_2));
title('Direct quadrature (16 points), integral 2')

num_levels = 8;
levels = linspace(-15, 1, num_levels+1); 
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
contour(real(Z), imag(Z), log10(abs(est2)), levels(2:end), 'k')
plot(complex(g(t16)), '.-r')
axis equal

%-------------------------------------------------------------------------
figure(3)
clf
pcolor(real(Z), imag(Z), log10(E16_3));
title('Direct quadrature (16 points), integral 3')

num_levels = 8;
levels = linspace(-15, 1, num_levels+1); 
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
contour(real(Z), imag(Z), log10(abs(est3)), levels(2:end), 'k')
plot(complex(g(t16)), '.-r')
axis equal

%-------------------------------------------------------------------------
figure(4)
clf
pcolor(real(Z), imag(Z), log10(E16_tot));
title('Direct quadrature (16 points), integral tot')

num_levels = 8;
levels = linspace(-15, 1, num_levels+1); 
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
contour(real(Z), imag(Z), log10(abs(esttot)), levels(2:end), 'k')
plot(complex(g(t16)), '.-r')
axis equal

disp('End!')

end

% ====================== MORE FUNCTIONS

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
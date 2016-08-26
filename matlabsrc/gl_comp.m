% Compare error for Laplace on starfish with trap. rule, GL and specialized
% GL

clear all; close all; clc

% ------------------- Set up parametrization ----------------------------
tau = @(t) (1+0.3*cos(5*t)).*exp(1i*t);
tau_prim = @(t) (-1.5*sin(5*t)+1i*(1+0.3*cos(5*t))).*exp(1i*t);
tau_primprim = @(t) exp(1i*t).*(-1-7.8*cos(5*t)-(3i)*sin(5*t));

zp = 1.5 + 1.5i;
RHS = @(x) imag(x.^2./(x-zp));

% ------------------- Fill domain with computation points ---------------
Nr_fill = 100;
r = linspace(0,0.8,Nr_fill)';
rr = linspace(0.8,1,round(Nr_fill/2))';
r = [r; rr];

Nt_fill = 64;
t = linspace(0,2*pi,Nt_fill)';
[R,T] = meshgrid(r,t);

z_starfish = R(:).*tau(T(:));
uk = RHS(z_starfish); %compute known exact solution

Z = R.*tau(T); 
%%

% ---------------- Compute using trapezoidal rule ------------------------
N = 256; %nbr discretization points on s
s = linspace(0,2*pi,N)';
z_trap = tau(s);
ub_trap = RHS(z_trap);

u_trap = laplace_trap(z_trap,tau_prim(s),tau_primprim(s),zp,z_starfish);
e_trap = abs(u_trap-uk);

E_trap = reshape(log10(e_trap),size(R));
figure()
pcolor(real(Z),imag(Z),E_trap)
shading interp
colormap parula
colorbar
title('10log error trap. rule')
caxis([-16 13])

%%
% ---------------- Compute using normal Gauss-Legendre quadrature -------
N_panels = 30; %nbr GL panels on starfish (16 points each)
panels = linspace(0,2*pi,N_panels+1);
z_gl = zeros(N_panels*16,1);
s_gl = z_gl;
w_gl = zeros(N_panels*16,1);
for i=1:N_panels
    [nodes,weights] = gaussleg(16,[panels(i) panels(i+1)]);
    s_gl((i-1)*16+1:16*i) = nodes;
    z_gl((i-1)*16+1:16*i) = tau(nodes);
    w_gl((i-1)*16+1:16*i) = weights;
end

[u_gl,mu_gl] = laplace_gl(z_gl,tau_prim(s_gl),tau_primprim(s_gl),zp,z_starfish,w_gl);
e_gl = abs(u_gl-uk);

E_gl = reshape(log10(e_gl),size(R));
figure()
pcolor(real(Z),imag(Z),E_gl)
shading interp
colormap parula
colorbar
title('10log error normal GL')
caxis([-16 13])

%% Calculate velocity in the domain using mu from normal and GL quadrature 

zDrop = z_gl;
WDrop = w_gl;
peDrop = panels;

% [u] = mex_spec_quad(zDrop,zpDrop,WDrop,peDrop,mu_gl)
[u] = mex_spec_quad(zDrop,WDrop,peDrop,mu_gl);

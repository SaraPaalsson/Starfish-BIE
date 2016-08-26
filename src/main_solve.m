% Solve Laplace's equation on starfish domain with Dirichlet b.c.

% close all; clear all; clc

% ----------------- Construct starfish domain ----------------------------
% Create starfish parametrization
% Obtan GL-16 nodes and weights
% Fill domain with computation points
s = 2/pi;%4/pi; %panel length = 0.2094
pl = 0.2094;
tau = @(t) (1+0.3*cos(5*(t+s))).*exp(1i*(t+s)); %starfish parametrization
taup = @(t) (-1.5*sin(5*(t+s))+1i*(1+0.3*cos(5*(t+s)))).*exp(1i*(t+s));
taupp = @(t) exp(1i*(t+s)).*(-1-7.8*cos(5*(t+s))-(3i)*sin(5*(t+s)));
% tau = @(t) cos(t) + 1i*sin(t); %circle parametrization
% taup = @(t) -sin(t) + 1i*cos(t);
% taupp = @(t) -cos(t) -1i*sin(t);

load 'glW.mat' %read in GL 16 and 32 weights

Npanels = 35; %nbr of panels on interface
N = Npanels*16; %nbr of discretization points on interface

[zDrops, wDrops,tpar,panels] = gl16(Npanels,tau);


%Test different resolutions here
[z,zplot,Rplot,Tplot] = fill_starfish(tau);
% K = 1; levelsK = 98; Nr_fill = 100;
% [z,relp] = fill_starfish2(zDrops,K, levelsK,Nr_fill);


%%
% ----------------- Set up problem ---------------------------------------
% Define boundary condition, Laplace's eq.
% Calculate known solution over the domain
zp = 1.5 + 1.5i;
RHS = @(x) imag(x.^2./(x-zp));

u_known = RHS(z);

% ----------------- Calculate density ------------------------------------
% Solve BIE to obtain density mu

mu_lapl = mubie_lapl(N,zDrops,taup(tpar),taupp(tpar),wDrops,RHS);

%%
% ----------------- Calculate u ------------------------------------------
% Compute u over the domain
% Use 16-GL when possible% As done in Trefethen
% Use special quadrature for points too close to the boundary

load 'IP1632.mat'

u = compu_lapl(N, mu_lapl, z, zDrops, taup(tpar), wDrops);
[uspec,z16,z32,zinterp] = specquad_lapl(u, mu_lapl, Npanels, panels, zDrops, ...
    taup(tpar), wDrops, z, tau, IP1, IP2, W16, W32,u_known);


%%

% ----------------- Plot solution and error ------------------------------
% plot with meshgrid variables

% rel_plot();
% rel_plot_av();
filled_plot();

%%
% % figure()
% % plot(zDrops,'.-')
% % hold on
% % plot(z(1:N),'r.')
% % plot(z(N+1:2*N),'g.'
% av_er = zeros(levelsK,1);
% for j=1:length(relp)
%     indie = ((j-1)*N+1:j*N);
%     era = abs(u_known(indie)-uspec(indie));
%    av_er(j) =  sum(era)/N;
% end
% %%
% 
% 
% 
% figure()
% plot(relp(1:end-1),av_er,'.-')
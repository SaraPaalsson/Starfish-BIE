function [z,zplot,Rplot,Tplot] = fill_starfish(tau)
% Fills the domain as described by tau with computation points. We
% discretize the radius r and the parametrization t. Create R and T with
% meshgrid for plotting later on
Nr_fill = 100;
r = linspace(0,0.9999,Nr_fill)';
% % rr = linspace(0.8,0.9999,Nr_fill)';
% % rr = rr(1:end-10);
% % r = [r; rr];

% % Nr_fill = 200;
% % r = logspace(0,0.9999,Nr_fill)';
% % r2 = (10-r)/10;


% Nr_fill = 6;
% relp = logspace(-5,0,Nr_fill)'; relp=relp(1:end-1);
% relp = flip(relp);
% r = 1-relp;


% Nt_fill = 128;
Nt_fill = 128;
t = linspace(0,2*pi,Nt_fill)'; %t=t(1:end-1);
[Rplot,Tplot] = meshgrid(r,t);

z = Rplot(:).*tau(Tplot(:));
% z = Rplot.*tau(Tplot);
% z = z(:);
zplot = Rplot.*tau(Tplot); 
end
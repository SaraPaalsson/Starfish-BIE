function [z,relp] = fill_starfish2(zDrops,K,levelsK,Nr_fill)
% [z,zplot,Rplot,Tplot]
% Nt_fill = length(tpar);
% t = tpar;


% Nr_fill = 4;
relp = logspace(-15,0,Nr_fill)'; relp = relp(1:end-1);
relp = flip(relp);
r = 1-relp;

% z = (1-relp)*zDrops(K);


z = r*zDrops';
z = z';
z = z(:);

% Nt_fill = 128;
% t = linspace(0,2*pi,Nt_fill)';
% [Rplot,Tplot] = meshgrid(r,t);


% 
% 
% z = (1-relp)*zDrops(K);

% K2 = K:1:K+levelsK;
% z = [];
% for i=1:length(K2)
%     z = [z; (1-relp(end-K2(i)))*zDrops];
% end



% Nr_fill = 10;
% rel_prox = logspace(-5,0,Nr_fill)';
% rel_proxflip = flip(rel_prox);
% 
% 
% b = tau(tpar(1));
% r = b - (rel_proxflip)*b;
% 
% z = r.*tau(tpar(1));

% Nr_fill = 100;
% r = linspace(0,0.9999,Nr_fill)';
% rr = linspace(0.8,0.9999,Nr_fill)';
% rr = rr(1:end-10);
% r = [r; rr];

% Nr_fill = 200;
% r = logspace(0,0.9999,Nr_fill)';
% r2 = (10-r)/10;

% r = zeros(Nt_fill,1);
% for j=1:Nt_fill
%     r((j-1)*Nr_fill+1) = tau(tpar(j)) + 10.^(rel_prox
% end
% 
% r = tau(tpar) + 10.^(rel_prox).*tau(tpar);
% 
% Nt_fill = 128;
% t = linspace(0,2*pi,Nt_fill)';
% [Rplot,Tplot] = meshgrid(r,t);
% 
% z = Rplot(:).*tau(Tplot(:));
% zplot = Rplot.*tau(Tplot); 
end
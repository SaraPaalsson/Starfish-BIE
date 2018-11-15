function errorest_mapping()
close all;
clear;
clc

defaultplotscript

% Interface we are regarding
% zfunc = @(ti) ti + 1i*cos(ti)% + 1i*sin(ti);
zfunc = @(ti) ti + 1i*ti.^2;
% zfunc = @(ti) ti; 
% panelp = [pi/4 1.2*pi/4];
panelp = [-1 1];
% panelp = [0 1];

% True gamma
g = @(ti) ones(size(ti));
% g = @(ti) sin(2*pi*ti);

% Discretise with GL, order n
n = 16;
[tn,wn] = gaussleg(n,[-1 1]);
tn = (tn+1)*(panelp(2)-panelp(1))/2 + panelp(1);
zinterf = zfunc(tn);  % Interface disc. points
wn = (panelp(2)-panelp(1))/2*wn;

% gn = g(zinterf);

figure(1); plot(real(zinterf),imag(zinterf),'.-','MarkerSize',20); hold on; grid on
title('Interface')

% Scale interface such that it has endpoints at -1 and 1
mid = (zfunc(panelp(2))+zfunc(panelp(1)))/2;
len = zfunc(panelp(2)) - zfunc(panelp(1));
zinterf_sc = 2*(zinterf-mid)/len;

gn = g(zinterf);

figure(1); plot(real(zinterf_sc),imag(zinterf_sc),'.-','MarkerSize',20)

% Define point(s) where to compute error
z0 = 0.25 + 1.25i;
figure(1); plot(real(z0),imag(z0),'*');

% Scale target point(s)
nz = 2*(z0-mid)/len;
figure(1); plot(real(nz),imag(z0),'o')

% Find coefficients for gamma:
ghat_l = comp_ghatl(n,tn,wn,gn);

% Test that the interpolant and ghat_l are correct by interpolating back at
% tn:
Pn_tn = zeros(n,1); lvec = (0:n-1)';
for j=1:n
    Pn_tn(j) = sum(ghat_l.*legendreP(lvec,zinterf_sc(j)));
end
Pn_tn = Pn_tn;%/((panelp(2)-panelp(1))/2);

figure(2); plot(tn,gn); hold on
plot(tn,Pn_tn,'*'); 
title('Gamma and interpolation (*)'); grid on

figure(3); semilogy(abs(Pn_tn-gn)+eps,'.-')
title('Interpolation error'); grid on

end

function ghat_l = comp_ghatl(n,ti,wi,gi)

ghat_l = zeros(n,1);
for l=0:n-1
    ghat_l(l+1) = (2*l+1)/2*sum(legendreP(l,ti).*wi.*gi);
end

end




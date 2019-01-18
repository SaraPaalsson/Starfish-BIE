%Test specquad stokes.

close all; clear; clc

Npanels = 1;

% Panel straight line, real
tau = @(t) (-pi + t)/pi+1i*(-pi + t)/pi;
taup = @(t) (1/pi+1i/pi)*ones(size(t));

% Obtan GL-16 nodes and weights for panel 
[zDrops, wDrops,tpar,panels] = gl16(Npanels,tau);
zpDrops = taup(tpar);
tau_panels = tau(panels);

z = 0.5i;

figure(1); plot(real(zDrops),imag(zDrops),'.-','MarkerSize',15); hold on;
axis equal; axis square 
plot(real(tau_panels),imag(tau_panels),'o')
figure(1); plot(real(z),imag(z),'k*')
%%
f = @(x) ones(size(x));

K = f(zDrops).*zpDrops./(zDrops-z);

Kf = @(t) f(tau(t)).*taup(t)./(tau(t)-z);

% figure(2); plot(zDrops,real(K),'.-'); hold on; plot(zDrops,imag(K),'.--')

IM = integral(Kf,0,2*pi);
IGL = sum(wDrops.*K);
% Icorr = log(1-z) - log(-1-z);
% errGL = abs(Icorr-IGL)
% errM = abs(Icorr-IM)
errGL_M = abs(IM-IGL)

mid2 = (tau_panels(2:end)+tau_panels(1:end-1))/2; %Midpoint panel
len2 = tau_panels(2:end)-tau_panels(1:end-1); %Length panel

tz = zDrops;
nzpan = 2*(tz-mid2)/len2;  %map all nodes on the panel to [-1,1]
nz = 2*(z-mid2)/len2; %map zk to nz (z0 in paper)

lg1 = log(1-nz);
lg2 = log(-1-nz); 
p32 = zeros(32,1);
p32(1) = lg1-lg2;
% div1 = 1/(1+nz);
% div2 = 1/(1-nz);
% q32 = zeros(32,1);
% q32(1) = -div1 - div2;

tzp = zpDrops;
tW = wDrops;            

load 'glW.mat' %read in GL 16 and 32 weights
load 'IP1632.mat'
% Interpolate to 32-point GL quadrature
tz32 = IPmultR(tz,IP1,IP2);
tzp32 = IPmultR(tzp,IP1,IP2);
nzpan32 = IPmultR(nzpan,IP1,IP2);

signc = -1;
for j=1:31 %Calculate pk:s...
    p32(j+1) = nz*p32(j) + (1-signc)/(j);
    signc = -signc;
%     q32(j+1) = nz*q32(j) + p32(j);
end
                    
p32coeff = vandernewton(nzpan32,p32,32);
% q32coeff = vandernewton(nzpan32,q32,32);

ISQ = sum(f(tz32).*p32coeff);
% errSP = abs(Icorr-ISQ)
errSP_M = abs(IM-ISQ)
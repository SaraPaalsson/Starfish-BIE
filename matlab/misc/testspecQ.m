%Test specquad stokes.

close all; clear; clc

Npanels = 1;

% % Panel straight line, real
% tau = @(t) (-pi + t)/pi;
% taup = @(t) 1/pi*ones(size(t));

% % Panel line, -1-i -> 1+i
% tau = @(t) (-pi + t)/pi + 1i*(-pi+t)/pi;
% taup = @(t) (1/pi + 1i/pi)*ones(size(t));

% Panel sin, 
tau = @(t) (-pi + t) + 1i*sin((-pi+t));%-1i*t;
taup = @(t) ones(size(t)) + 1i*cos(-pi+t);


% Obtan GL-16 nodes and weights for panel 
[zDrops, wDrops,tpar,panels] = gl16(Npanels,tau);
zpDrops = taup(tpar);
tau_panels = tau(panels);

z = .7i;

figure(1); plot(real(zDrops),imag(zDrops),'.-','MarkerSize',15); hold on;
axis equal; axis square 
plot(real(tau_panels),imag(tau_panels),'o')
figure(1); plot(real(z),imag(z),'k*')


f = @(x) x-z; %x-0.5*z; %ones(size(x)); %sin((x));
fDrops = f(zDrops);

K = f(zDrops).*zpDrops./(zDrops-z).^2;

Kf = @(t) f(tau(t)).*taup(t)./(tau(t)-z).^2;


% figure(); plot(real(zDrops),real(fDrops),'.-'); hold on; plot(real(zDrops),imag(fDrops),'.--')

IM = integral(Kf,0,2*pi);
IGL = sum(wDrops.*K);
% Icorr_1 = 1/(z-1)-1/(z+1);
% Icorr_2 = (0.5 - 1i*0.5)*log(-nz+(1+1i)*1)-(0.5 - 1i*0.5)*log(-nz+(1+1i)*(-1));
% Icorr = Icorr_2;
% errGL = abs(Icorr-IGL)
% errM = abs(Icorr-IM)
errGL_M = abs(IM-IGL)

mid2 = (tau_panels(2:end)+tau_panels(1:end-1))/2; %Midpoint panel
len2 = tau_panels(2:end)-tau_panels(1:end-1); %Length panel
nz = 2*(z-mid2)/len2; %map zk to nz (z0 in paper)

tz = zDrops;
nzpan = 2*(tz-mid2)/len2;  %map all nodes on the panel to [-1,1]

figure(2); plot(real(nzpan),imag(nzpan),'.-','MarkerSize',15); hold on;
axis equal; axis square 
plot([-1 1],[0 0],'o')
plot(real(nz),imag(nz),'k*')
p32 = zeros(32,1);
q32 = zeros(32,1);

lg1 = log(1-nz);
lg2 = log(-1-nz); 
div1 = 1/(1+nz);
div2 = 1/(1-nz);
q32(1) = -div1 - div2;

% lg1 = lg1 - pi*1i; % Above real line, under panel
% lg2 = lg2 + pi*1i; % Above real line, under panel

% lg1 = lg1 + pi*1i; % Below real line, above panel
% lg2 = lg2 - pi*1i; % Below real line, above panel

p32(1) = lg1-lg2;


tzp = zpDrops;
tW = wDrops;            

load 'glW.mat' %read in GL 16 and 32 weights
load 'IP1632.mat'
% Interpolate to 32-point GL quadrature
tz32 = IPmultR(tz,IP1,IP2);
tzp32 = IPmultR(tzp,IP1,IP2);
nzpan32 = IPmultR(nzpan,IP1,IP2);

alpha = 2/len2;
% alpha = len2/2;
 
signc = -1;
for j=1:31 %Calculate pk:s...
    p32(j+1) = nz*p32(j) + (1-signc)/(j);
    signc = -signc;
    q32(j+1) = nz*q32(j) + p32(j);
end
                    
p32coeff = vandernewton(nzpan32,p32,32);
q32coeff = vandernewton(nzpan32,q32,32);

ISQ = alpha*sum(f(tz32).*q32coeff);
errSP_M = abs(IM-ISQ)

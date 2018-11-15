%Test specquad stokes.

close all; clear; clc

Npanels = 1;
% 
% % Panel straight line, real
% tau = @(t) (-pi + t)/pi;
% taup = @(t) 1/pi*ones(size(t));

% Panel line, -1-i -> 1+i
tau = @(t) (-pi + t)/pi-0.1i + 1i*(-pi+t)/pi;
taup = @(t) (1/pi + 1i/pi)*ones(size(t));

% % Panel sin, 
% tau = @(t) (-0.1 + t) + 1i*sin((-0.1+t));
% taup = @(t) ones(size(t)) + 1i*cos(-0.1+t);


% Obtan GL-16 nodes and weights for panel 
[zDrops, wDrops,tpar,panels] = gl16(Npanels,tau);
zpDrops = taup(tpar);
nDrops = -1i*zpDrops./abs(zpDrops);
tau_panels = tau(panels);

z = 1i;

figure(1); plot(real(zDrops),imag(zDrops),'.-','MarkerSize',15); hold on;
axis equal; axis square 
plot(real(tau_panels),imag(tau_panels),'o')
figure(1); plot(real(z),imag(z),'k*')

% f = @(x) sin(x); %ones(size(x)); %sin((x));
% fDrops = f(zDrops);
% K = f(zDrops).*zpDrops./(zDrops-z).^2;
% Kf = @(t) f(tau(t)).*taup(t)./(tau(t)-z).^2;

K = conj(sin(zDrops)).*imag(conj(zDrops-z).*zpDrops)./(conj(zDrops-z).^2);
Kf = @(t) conj(sin(tau(t))).*imag(conj(tau(t)-z).*taup(t))./(conj(tau(t)-z).^2);


% figure(); plot(real(zDrops),real(fDrops),'.-'); hold on; plot(real(zDrops),imag(fDrops),'.--')

IM = (integral(Kf,0,2*pi));
IGL = (sum(wDrops.*K));
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
nz32 = -1i*tzp32./abs(tzp32);
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

f1 = sin(tz32).*conj(nz32).^2;
f2 = sin(tz32).*conj(tz32-z);

I1 = sum(f1.*p32coeff);
I2 = alpha*sum(f2.*q32coeff);


ISQ = 1i/2*conj(I1) + 1i/2*conj(I2);
errSP_M = abs(IM-ISQ)

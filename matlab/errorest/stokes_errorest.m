close all;
clear; 
clc

load('IP1632.mat')
load('glW.mat')

% Do error estimates for Stokes! (using Ludvigs paper)

% For tests, define a nice density function
mu = @(t) ones(size(t));

% Also for tests, do integrals on flat Gamma on [-1,1]
n = 16;
[zn,wn] = gaussleg(n,[-1 1]);
znp = ones(size(zn));
mun = mu(zn);
nt = -1i*znp./abs(znp); %

% Evaluation point in domain
a = 0; bv = linspace(0.001,1,100)'; %OBS also works for a != 0
% bv = 1;
z0 = a + bv*1i;

% Define approximation of kn
kn = @(x) 2*pi./(abs(abs(x) + abs(sqrt(x.^2-1))).^(2*n+1));
tcirc = linspace(0,2*pi,1001)'; tcirc = tcirc(1:end-1);
zcirc = @(r) r*cos(tcirc) + 1i*r*sin(tcirc);
zpcirc = @(r) -r*sin(tcirc) + 1i*r*cos(tcirc);
wcirc = tcirc(2)-tcirc(1);

% Compute integral with quadrature
I1 = zeros(length(bv),1); I2 = I1; mu2max = I2;
est1 = I1; err1 = I1; I1_spec = I1; err1_spec = I1;
I2_spec = I2; err2_spec = I2; est2 = est1;
I3 = I2; est3 = est1; I3_spec = I3; err3_spec = err2_spec;
I4 = I2; I4_spec = I2_spec; est4=est1; err4_spec = err2_spec;
for j=1:length(bv)
    I1(j) = sum(mun.*wn.*imag(1./(zn-z0(j))).*znp);
    I2(j) = conj(sum(mun.*conj(nt).^2.*wn./(zn-z0(j)).*znp));
    I3(j) = conj(sum(mun.*wn.*conj(zn-z0(j)).*znp./((zn-z0(j)).^2)));
    I4(j) = sum(mun.*wn.*imag(znp.*conj(zn-z0(j)))./(conj(zn-z0(j)).^2));
    
    % Compute exact solution for integral of f(z), f(t) = 1/(t-z0)
    b = imag(z0(j));
    I1r = 0.5*log(b.^2+(1-a)^2) - 0.5*log(b.^2+(-1-a)^2);
    I1i = atan((1-a)./b) - atan((-1-a)./b);
    I1_correct = I1i;

    % Compute solution with special quadrature,
    p32 = zeros(32,1);
    p32(1) = log(1-z0(j))-log(-1-z0(j));
    q32 = zeros(32,1);
    q32(1) = -1/(1+z0(j)) - 1/(1-z0(j));
    tmu32 = IPmultR(mun,IP1,IP2);
    tz32 = IPmultR(zn,IP1,IP2);
    tzp32 = IPmultR(znp,IP1,IP2);
    
    gamma = 2/2;
    signc = -1;
    for jk=1:31 %Calculate pk:s...
        p32(jk+1) = z0(j)*p32(jk) + (1-signc)/(jk);
        signc = -signc;
        q32(jk+1) = z0(j)*q32(jk) + p32(jk);
    end
                    
    ntz = -1i*tzp32./abs(tzp32); %OK
    f1 = tmu32.*conj(ntz).^2; %OK
    f2 = tmu32.*conj(tz32-z0(j)); %OK
    p32coeff = vandernewton(tz32,p32,32);
    q32coeff = vandernewton(tz32,q32,32);
                                        
    new1 = (transpose(tmu32)*imag(p32coeff(1:end)));
    new2 = conj(sum(f1.*p32coeff(1:end)));
    new3 = conj(gamma*sum(f2.*q32coeff(1:end)));

    I1_spec(j) = new1;
    err1_spec(j) = abs(I1(j)-new1);
    I2_spec(j) = new2;
    err2_spec(j) = abs(I2(j)-new2);
    I3_spec(j) = new3;
    err3_spec(j) = abs(I3(j)-new3);
    I4_spec(j) = 1i/2*new2 + 1i/2*new3;
    err4_spec(j) = abs(I4(j)-I4_spec(j));
    
    % Compute estimate
    mumax1 = max(abs(mun.*znp));
    est1(j) = abs(mumax1*2*pi./(abs(z0(j) + sqrt(z0(j).^2-1)).^(2*n+1)));
    mumax2 = max(abs(mun.*conj(nt).^2).*znp);
    est2(j) = abs(mumax2*2*pi./(abs(z0(j) + sqrt(z0(j).^2-1)).^(2*n+1)));
    mumax3 = max(abs(mun.*znp));
    est3(j) = abs(pi./((z0(j) + sqrt(z0(j)^2-1))^(2*n+1))...
        *((2*n+1)^2/(sqrt(z0(j)^2-1)^2) + z0(j)*(2*n+1)/(sqrt(z0(j)^2-1)^3)));
    
%     R = max(abs(z0(j)-1),abs(z0(j)-(-1)));
%     zcircle = zcirc(R);
%     zpcircle = zpcirc(R);
%     
%     f3 = kn(zcircle).*zpcircle./((zcircle-z0(j)).^3);
%     intR3 = abs(1/(2i*pi)*sum(f3.*wcirc));
%     est3(j) = abs(intR3 - est3(j));
%     

    err1(j) = abs(I1(j)-I1_correct);
end

figure(1);clf;
subplot(131);
plot(bv,real(I1),'b-','DisplayName','I1 comp RE'); hold on;
plot(bv,real(I1_spec),'r-','DisplayName','I1 spec. IM');
plot(bv,imag(I1),'b--','DisplayName','I1 comp RE'); hold on;
plot(bv,imag(I1_spec),'r--','DisplayName','I1 spec. IM');
h=legend('toggle'); set(h,'Location','Southwest')
title('I1')
subplot(132)
plot(bv,real(I2),'b-','DisplayName','I2 comp RE'); hold on;
plot(bv,real(I2_spec),'r-','DisplayName','I2 spec. IM');
plot(bv,imag(I2),'b--','DisplayName','I2 comp RE'); hold on;
plot(bv,imag(I2_spec),'r--','DisplayName','I2 spec. IM');
h=legend('toggle'); set(h,'Location','Southwest')
title('I2')
subplot(133)
plot(bv,real(I3),'b-','DisplayName','I3 comp RE'); hold on;
plot(bv,real(I3_spec),'r-.','DisplayName','I3 spec. IM');
plot(bv,imag(I3),'b--','DisplayName','I3 comp RE'); hold on;
plot(bv,imag(I3_spec),'r--','DisplayName','I3 spec. IM');
h=legend('toggle'); set(h,'Location','Southwest')
title('I3')

figure(2); clf;
subplot(131);
% semilogy(bv,err1,'.-','DisplayName','Rel.err'); hold on; grid on
semilogy(bv,err1_spec,'DisplayName','Rel.err to spec'); hold on; grid on
semilogy(bv,est1,'k--','DisplayName','Est.')
xlabel('$\Im(z_0)$')
h=legend('toggle'); set(h,'Location','Southwest')
title('I1')
subplot(132);
semilogy(bv,err2_spec,'DisplayName','Rel.err to spec'); hold on; grid on
semilogy(bv,est2,'k--','DisplayName','Est.')
xlabel('$\Im(z_0)$')
h=legend('toggle'); set(h,'Location','Southwest')
title('I2')
subplot(133);
semilogy(bv,err3_spec,'DisplayName','Rel.err to spec'); hold on; grid on
semilogy(bv,est3,'k--','DisplayName','Est.')
xlabel('$\Im(z_0)$')
h=legend('toggle'); set(h,'Location','Southwest')
title('I3')

figure(3); clf;
semilogy(bv,abs(err3_spec-est3),'k:')
grid on



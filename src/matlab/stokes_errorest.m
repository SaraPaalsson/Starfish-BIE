close all;
clear; 
clc

% Do error estimates for Stokes! (using Ludvigs paper)

% For tests, define a nice density function
mu = @(t) ones(size(t));

% Also for tests, do integrals on flat Gamma on [-1,1]
n = 16;
[zn,wn] = gaussleg(n,[-1 1]);
mun = mu(zn);

% Evaluation point in domain
a = 0; b = linspace(0.05,0.5,100)'; %OBS also works for a != 0
z0 = a + b*1i;

% figure(1); clf; 
% plot(zn,zeros(size(zn)),'k.-','MarkerSize',15); hold on
% plot(z0,'r*')
% title('Interface $\Gamma$'); grid on

% Compute integral with quadrature
I1 = zeros(length(b),1);
for j=1:length(b)
    I1(j) = sum(wn.*imag(1./(zn-z0(j))));
end

% Compute exact solution for integral of f(z), f(t) = 1/(t-z0)
I1r = 0.5*log(b.^2+(1-a)^2) - 0.5*log(b.^2+(-1-a)^2);
I1i = atan((1-a)./b) - atan((-1-a)./b);
% I1_correct = I1r + 1i*I1i;
I1_correct = I1i;


% Compute estimate
est = 2*pi./(abs(z0 + sqrt(z0.^2-1)).^(2*n+1));

err = abs(I1-I1_correct);

figure(2); clf;
plot(b,err,'.-','DisplayName','Rel.err'); hold on; grid on
plot(b,est,'k--','DisplayName','Est.')
xlabel('$\Im(z_0)$')
legend('toggle')
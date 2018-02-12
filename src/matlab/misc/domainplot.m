% Domain plot
% Make example plot of domain for single drop

defaultplotscript
addpath(genpath('~/Documents/MATLAB/matlab2tikz'))

s = 0;
a = 0.3;
tau = @(t) (1+a*cos(5*(t+s))).*exp(1i*(t+s)); %starfish parametrization
taup = @(t) -(1i*tau(t) + exp(1i*(t+s)).*(-5*a*sin(5*(t+s))));

N = 128;
alpha = linspace(0,2*pi,N+1)';

z = tau(alpha);
zp = taup(alpha);
zp = zp./abs(zp);

%% Plot interface tangent and normal
figure(1); clf
plot(z,'k','LineWidth',3);  
hold on; axis equal
ylim([-1.9 1.9])
set(gca,'Visible','Off')

% Plot tangential vector on omega 1 (z3)
ind = 21;
q = quiver(real(z(ind)),imag(z(ind)),real(zp(ind)),imag(zp(ind)));
q.Color = 'black';
q.LineWidth = 1.5;
q.MaxHeadSize = 1;
q.AutoScaleFactor = 0.7;
str = '$\mathbf{s}$';
t = text(0.7,1.25,str);
set(t,'interpreter','latex','FontSize',20)

% Plot normal vector
n = -1i*zp;
q = quiver(real(z(ind)),imag(z(ind)),real(n(ind)),imag(n(ind)));
q.Color = 'black';
q.LineWidth = 1.5;
q.MaxHeadSize = 1;
q.AutoScaleFactor = 0.7;
str = '$\mathbf{n}$';
t = text(0.1,0.95,str);
set(t,'interpreter','latex','FontSize',20)

str = '$\mathbf{\Omega_0}$';
t = text(1.5,-1.5,str);
set(t,'interpreter','latex','FontSize',20)


str = '$\mathbf{\Omega}$';
t = text(0,0,str);
set(t,'interpreter','latex','FontSize',20)


str = '$\mathbf{\Gamma}$';
t = text(-1.2,-0.35,str);
set(t,'interpreter','latex','FontSize',20)


% matlab2tikz

%% Plot interface with boundary or domain discretization
xl = -1.5; xu = -xl;
yl = -1.5; yu = -yl;
dt = 0.1;

Ind = 1:2:N;

figure(2); clf
hold on; axis equal
grid on
xlim([xl xu])
xticks(xl:dt:xu)
ylim([yl yu])
yticks(yl:dt:yu)
% set(gca,'GridColor','k','LineWidth',0.5,'GridAlpha',0.5)
plot(z,'LineWidth',3); 

set(gca,'Visible','Off')
% plot(z(Ind),'k.','MarkerSize',20)

%% Plot composite discretization

Npanels = 10;
[z, ~,alpha,panels] = gl16(Npanels,tau);

zp = taup(alpha);

pz = tau(panels);
pz_p = taup(panels);
pz_n = 1i*pz_p./abs(pz_p);

k = 0.1;

figure(3); clf

atmp = linspace(0,2*pi,400)';
ztmp = tau(atmp);
plot(ztmp,'LineWidth',3)
hold on; axis equal
for j=1:length(pz)
    plot(complex([pz(j)-k*pz_n(j) pz(j)+k*pz_n(j)]),'k-')
end
set(gca,'Visible','Off')
xlim([-1.5 1.5])
ylim([-1.4 1.4])

figure(4); clf
plot(ztmp,'LineWidth',3)
hold on; axis equal
for j=1:length(panels)
    plot(complex([pz(j)-k*pz_n(j) pz(j)+k*pz_n(j)]),'k-','LineWidth',3)
end
% plot(z(1:16),'.','MarkerSize',35)
alphaequi = linspace(0,2*pi,Npanels*16+1)';
zequi = tau(alphaequi);
plot(zequi(1:16),'.','MarkerSize',35)
xlim([0.4 1.4])
ylim([-0.2 0.6])
% xlim([-1.4 1.4])
% ylim([-1.4 1.4])
set(gca,'Visible','Off')
publication_fig









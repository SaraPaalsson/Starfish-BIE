close all; clear all; clc

%% Solve BIE for starfish to compute complex density

% ----------------- Construct starfish domain ----------------------------
% Create starfish parametrization
% Obtan GL-16 nodes and weights
% Fill domain with computation points
s = 0;
tau = @(t) (1+0.3*cos(5*(t+s))).*exp(1i*(t+s)); %starfish parametrization
taup = @(t) (-1.5*sin(5*(t+s))+1i*(1+0.3*cos(5*(t+s)))).*exp(1i*(t+s));
taupp = @(t) exp(1i*(t+s)).*(-1-7.8*cos(5*(t+s))-(3i)*sin(5*(t+s)));


% load 'glW.mat' %read in GL 16 and 32 weights

%Npanels = 35; %nbr of panels on interface
 Npanels = 70;
N = Npanels*16; %nbr of discretization points on interface


[zDrops, wDrops,tpar,panels] = gl16(Npanels,tau);

% ----------------- Set up problem ---------------------------------------
% Define boundary condition, Laplace's eq.
%zp = 1.5 + 1.5i;
%RHS = @(x) imag(x.^2./(x-zp));

zsrc1 = 1.5+1.5i;
zsrc2 = -0.25+1.5i;
zsrc3 = -0.5-1.5i;
RHS = @(x) real( 1 ./ (x-zsrc1) + 1 ./ (x-zsrc2) + 1 ./ (x-zsrc3) );

% ----------------- Calculate density ------------------------------------
% Solve BIE to obtain density mu
mu_lapl = mubie_lapl(N,zDrops,taup(tpar),taupp(tpar),wDrops,RHS);


%% Fill starfish with points to compute quadrature errors

% Fill the domain with computational points, discretized in r and t
% Grid coarse far from and fine close to boundary
% Compute only on the first quadrant
%nbrR = 300;
%nbrT = 300;
% nbrR = 2000;
% nbrT = 2500;
nbrR = 500;
nbrT = 500;

%R1 = 0.5; %Where to go from coarse to fine grid
R1 = 0.4;

r = [linspace(0,R1,10) linspace(R1, 0.999,nbrR)]';
% r = [linspace(R1, 0.999,nbrR)]';
t = linspace(0,pi/2,nbrT)';

[Rplot,Tplot] = meshgrid(r,t);
z = Rplot(:).*tau(Tplot(:));
zplot = Rplot.*tau(Tplot);

% Calculate known solution over the domain
u_known = RHS(z);


%% Compute u with normal and special quadrature
load 'IP1632.mat'
load 'glW.mat' %read in GL 16 and 32 weights

% ----------------- Calculate u ------------------------------------------
% Compute u over the domain
% Use 16-GL when possible
% Use special quadrature for points too close to the boundary
disp('Compute u normal quadrature')
tic
u = compu_lapl(N, mu_lapl, z, zDrops, taup(tpar), wDrops);
toc

% save('normquadu_high','u')

%%
% load('normquadu_high')

%disp('Compute u special quadrature')
%tic
%[uspec] = specquad_lapl(u, mu_lapl, Npanels, tau(panels), zDrops, ...
%    taup(tpar), wDrops, z, IP1, IP2, W16, W32);
%toc

 disp('Compute u special quadrature MEX')
 tic
 [uspec,~] = mex_saraspecquad(u, mu_lapl, tau(panels), zDrops, taup(tpar), wDrops, z);
 toc
% %
% save('specquadu_high','uspec')


%----------------------------------------------------
%% Compute error estimates

% disp('Compute estimates')
% 
% % Compute error estimate on grid
% errest = error_estL(z,tau(panels),zDrops,Npanels,mu_lapl,wDrops,taup(tpar));
% 

%%
% load('normquadu_high')
% load('specquadu_high')

% -----------------------------------------------------------
%% Compute errors and plot
close all;

disp('Compute  errors')
%filled_plot()

relnorm = norm(u_known,Inf);
er = abs(u_known-u)/relnorm;
Er = reshape(log10(er),size(Rplot));
er_spec = abs(u_known-uspec)/relnorm;
Erspec = reshape(log10(er_spec),size(Rplot));

% errest2 = {errest/relnorm reshape(errest/relnorm,size(Rplot))};

error = {reshape(er,size(Rplot))  reshape(er_spec,size(Rplot))};

% save('errquad_high','error','errest2')

%save('../HIGHRESdata35','u','uspec','error','errest2','zplot');
%save('../HIGHRESdata70','u','uspec','error','errest2');
%save('../LOWRESdata70','u','uspec','error','errest2');
%save('../LOWRESdata35','u','uspec','error','errest2');

%%

disp('Plot!')

%spec_str1 = sprintf('Quadrature error, %d panels, normal quadrature',Npanels);
%spec_str2 = sprintf('Quadrature error, %d panels, special quadrature',Npanels);
spec_str1 = '10log error normal quad., 35 panels'; spec_str2 = '10log error special quad., 35 panels';
titstr = {spec_str1 spec_str2};


% % For contour plot
levels = -15:3:-3; %????
%zoombox = [0.5 1.1 0.2 0.43]; %For 70 panels 
zoombox = [0.45 1.1 0.1 0.43];
%zoombox = [0.48 1.1 0.15 0.45];
drawbox = @(x,y) plot([x(1) x(2) x(2) x(1) x(1)],[y(1) y(1) y(2) y(2) y(1)],'-k');
tplot = linspace(0,2*pi,1000);


for i=1:2
    sfigure(i);
    clf
    publication_fig
    pcolor(real(zplot),imag(zplot),log10(error{i}))
    shading flat
    % shading interp
    colormap parula
    colorbar
    caxis([-18 -2])
    %ylabel('10log error normal G.-L.')
    set(gca,'yaxislocation','right');
    set(gca,'YTicklabel',[])
    set(gca,'XTicklabel',[])
    %    title(titstr{i})
    hold on
    plot(real(tau(tplot)), imag(tau(tplot)),'k')
    plot(real(zDrops),imag(zDrops),'.k','MarkerSize',3)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis([0 1.3 0 1.3])
    caxis([-15 0])
%     drawbox(zoombox(1:2),zoombox(3:4));
    publication_fig
    box on
end



sfigure(3);
clf
publication_fig
pcolor(real(zplot),imag(zplot),log10(error{2}))
shading flat
% shading interp
colormap parula
colorbar
caxis([-18 -2])
%ylabel('10log error normal G.-L.')
set(gca,'yaxislocation','right');
set(gca,'YTicklabel',[])
set(gca,'XTicklabel',[])
%    title(titstr{i})
hold on
plot(real(tau(tplot)), imag(tau(tplot)),'k')
% plot(real(zDrops),imag(zDrops),'.k','MarkerSize',4)
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
axis([0 1.3 0 1.3])
caxis([-15 0])
%    drawbox(zoombox(1:2),zoombox(3:4));
publication_fig
box on




%titstr = {'Level contours of log_{10} e(z), normal quadrature' 'Level contours of log_{10} e(z), special quadrature'};
titstr = {'Level curves, normal quad.' 'Level curves, special quad.'};
for i=4:5
    sfigure(i);
    clf;
    publication_fig
    plot(real(tau(tplot)), imag(tau(tplot)),'k')
    hold on
    contour(real(zplot),imag(zplot),log10(error{1}),levels,'k')
%     if i==5
%         contour(real(zplot),imag(zplot),log10(abs(errest2{2})),levels,'r-')
%     end
    plot(real(zDrops),imag(zDrops),'.k')
    shading flat
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis(zoombox)
    axis([0 1.3 0 1.3])
    caxis([-15 0])
    axis(zoombox)
    %     title(titstr{i-2})
    publication_fig
    box on
end

% levels = [-12 -9 -6];
for i=6:6
    sfigure(i);
    clf
    publication_fig
    plot(real(tau(tplot)), imag(tau(tplot)),'k')
    hold on
    contour(real(zplot),imag(zplot),log10(error{2}),levels,'k')
    %     contour(real(zplot),imag(zplot),log10(abs(errest2{2})),levels,'r-')
    plot(real(zDrops),imag(zDrops),'.k')
    shading flat
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis(zoombox)
    axis([0 1.3 0 1.3])
    caxis([-15 0])
    axis(zoombox)
    %     title(titstr{i-2})
    publication_fig
    box on
end

%%
% sfigure(1)
% print -dpng -r400 '../graphics/starfish_error_normq'
% %crop('../graphics/starfish_error_normq');
% sfigure(2)
% print -dpng -r400 '../graphics/starfish_error_interpq'
% %crop('../graphics/starfish_error_interpq');
% sfigure(3)
% print -dpng -r400 '../graphics/starfish_contour_normq'
% %crop('../graphics/starfish_contour_normq');
% sfigure(4)
% print -dpng -r400 '../graphics/starfish_contour_interpq'
% %crop('../graphics/starfish_contour_interpq');

% disp('Save plots')
% 
% % name = {['../graphics/filled_error_panels' num2str(Npanels)], ...
% %         ['../graphics/fillederror_SQbox_panels' num2str(Npanels)], ...
% %         ['../graphics/fillederror_SQ_panels' num2str(Npanels)], ...
% %         ['../graphics/contour_panels' num2str(Npanels)], ...
% %         ['../graphics/contour__LC_panels' num2str(Npanels)], ...
% %         ['../graphics/contour_SQ_panels' num2str(Npanels)]};
% 
% 
% 
% name = {['filled_error_panels' num2str(Npanels)], ...
%         ['fillederror_SQbox_panels' num2str(Npanels)], ...
%         ['fillederror_SQ_panels' num2str(Npanels)], ...
%         ['contour_panels' num2str(Npanels)], ...
%         ['contour__LC_panels' num2str(Npanels)], ...
%         ['contour_SQ_panels' num2str(Npanels)]};
% 
% for i=1:6
% %     sfigure(1)
% %     tmp = char(name{i})
% %     print -dpng -r400 tmp
%     print(i,name{i},'-dpng','-r400')
% %	crop(name{i});
% end

disp('Done!')

% name1 = ['../graphics/fillederror_panels' num2str(Npanels)]
% print -dpng -r400 name1
% printfig(1,name1)

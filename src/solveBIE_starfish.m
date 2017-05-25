%% Solve BIE for starfish to compute complex density
close all; clear all; clc

addpath('mex')

% ----------------- Set parameters ----------------------------
res_interf = 'low'; %low, high
res_domain = 'low'; %superlow, low, high
res = struct(); %result struct 
savePlots = 0; %if save plots
savedata = 0;

% ----------------- Setup domain  ------------------------------
[dom] = main_init(res_interf,res_domain);
res.dom = dom;

% ----------------- Set up problem -----------------------------
% Define boundary condition, Laplace's eq. --> Provides exact solution!
zsrc1 = 1.5+1.5i;
zsrc2 = -0.25+1.5i;
zsrc3 = -0.5-1.5i;
RHS = @(x) real( 1 ./ (x-zsrc1) + 1 ./ (x-zsrc2) + 1 ./ (x-zsrc3) );
res.RHS = RHS;


% ----------------- Calculate density ------------------------------------
% Solve BIE to obtain density mu
mu_lapl = mubie_lapl(dom.N,dom.zDrops,dom.taup(dom.tpar), ...
    dom.taupp(dom.tpar),dom.wDrops,RHS);
res.mu = mu_lapl;


% Calculate known solution over the domain
u_known = RHS(dom.z);
res.uknown = u_known;

% Compute u with normal and special quadrature
load 'IP1632.mat'
load 'glW.mat' %read in GL 16 and 32 weights

% ----------------- Calculate u ------------------------------------------
% Compute u over the domain
% Use 16-GL when possible
% Use special quadrature for points too close to the boundary
disp('Compute u normal quadrature')
tic
u = compu_lapl(dom.N, mu_lapl, dom.z, dom.zDrops, dom.taup(dom.tpar), dom.wDrops);
toc

if savedata 
    savestr = ['../results/normquadu_' res_domain];
    save(savestr,'u')
end
res.u = u;

% disp('Compute u special quadrature')
% tic
% [uspec] = specquad_lapl(u, mu_lapl, dom.Npanels, dom.tau(dom.panels), dom.zDrops, ...
%     dom.taup(dom.tpar), dom.wDrops, dom.z, IP1, IP2, W16, W32);
% toc
disp('Compute u special quadrature MEX')
tic
[uspec,~] = mex_saraspecquad(u, mu_lapl, dom.tau(dom.panels), dom.zDrops, dom.taup(dom.tpar), dom.wDrops, dom.z);
toc

if savedata
    savestr = ['../results/specquadu_' res_domain];
    save(savestr,'uspec')
end
res.uspec = uspec;

%----------------------------------------------------
% Compute error estimates

disp('Compute estimates')
% Compute error estimate on grid
errest = error_estL(dom.z,dom.tau(dom.panels),dom.zDrops,dom.Npanels, ...
    mu_lapl,dom.wDrops,dom.taup(dom.tpar)); 

savestr = ['../results/errorest_' res_domain];
save(savestr,'errest')
res.errest = errest;

%
% -------------------------
% 
% load('normquadu_high')
% load('specquadu_high')
% load('errorest_high')
% -------------------------

% Compute errors and plot
close all;

disp('Compute  errors')
%filled_plot()

relnorm = norm(u_known,Inf);
er = abs(u_known-u)/relnorm;
Er = reshape(log10(er),size(dom.zplot));
er_spec = abs(u_known-uspec)/relnorm;
Erspec = reshape(log10(er_spec),size(dom.zplot));

error = {reshape(er,size(dom.zplot))  reshape(er_spec,size(dom.zplot))};

if savedata
    savestr = ['../results/error_' res_domain];
    save(savestr,'error')
end
res.error = error;

%----------------------------------------------------
% Plot
disp('Plot!')
spec_str1 = '10log error normal quad., 35 panels'; spec_str2 = '10log error special quad., 35 panels';
titstr = {spec_str1 spec_str2};

% For contour plot
levels = -15:3:-3; 
zoombox = [0.45 1.1 0.1 0.43];
drawbox = @(x,y) plot([x(1) x(2) x(2) x(1) x(1)],[y(1) y(1) y(2) y(2) y(1)],'-k');
tplot = linspace(0,2*pi,1000);

for i=1:2
    sfigure(i);
    clf
    publication_fig
    pcolor(real(dom.zplot),imag(dom.zplot),log10(error{i}))
    shading flat
    % shading interp
    colormap parula
    colorbar
    caxis([-18 -2])
    set(gca,'yaxislocation','right');
    set(gca,'YTicklabel',[])
    set(gca,'XTicklabel',[])
    hold on
    plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
    plot(real(dom.zDrops),imag(dom.zDrops),'.k','MarkerSize',3)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis([0 1.3 0 1.3])
    caxis([-15 0])
    publication_fig
    box on
end

sfigure(3);
clf
publication_fig
pcolor(real(dom.zplot),imag(dom.zplot),log10(error{2}))
shading flat
% shading interp
colormap parula
colorbar
caxis([-18 -2])
set(gca,'yaxislocation','right');
set(gca,'YTicklabel',[])
set(gca,'XTicklabel',[])
hold on
plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
axis([0 1.3 0 1.3])
caxis([-15 0])
publication_fig
box on

titstr = {'Level curves, normal quad.' 'Level curves, special quad.'};
for i=4:5
    sfigure(i);
    clf;
    publication_fig
    plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
    hold on
    contour(real(dom.zplot),imag(dom.zplot),log10(error{1}),levels,'k')
    plot(real(dom.zDrops),imag(dom.zDrops),'.k')
    shading flat
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis(zoombox)
    axis([0 1.3 0 1.3])
    caxis([-15 0])
    axis(zoombox)
    publication_fig
    box on
end

for i=6:6
    sfigure(i);
    clf
    publication_fig
    plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
    hold on
    contour(real(dom.zplot),imag(dom.zplot),log10(error{2}),levels,'k')
    plot(real(dom.zDrops),imag(dom.zDrops),'.k')
    shading flat
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    axis(zoombox)
    axis([0 1.3 0 1.3])
    caxis([-15 0])
    axis(zoombox)
    publication_fig
    box on
end


%----------------------------------------------------
% Save plots if wanting to
if savePlots
    % disp('Save plots')
    sfigure(1)
    print -dpng -r400 '../graphics/starfish_error_normq'
    %crop('../graphics/starfish_error_normq');
    sfigure(2)
    print -dpng -r400 '../graphics/starfish_error_interpq'
    %crop('../graphics/starfish_error_interpq');
    sfigure(3)
    print -dpng -r400 '../graphics/starfish_contour_normq'
    %crop('../graphics/starfish_contour_normq');
    sfigure(4)
    print -dpng -r400 '../graphics/starfish_contour_interpq'
    %crop('../graphics/starfish_contour_interpq');
    
    name = {['filled_error_panels' num2str(Npanels)], ...
        ['fillederror_SQbox_panels' num2str(Npanels)], ...
        ['fillederror_SQ_panels' num2str(Npanels)], ...
        ['contour_panels' num2str(Npanels)], ...
        ['contour__LC_panels' num2str(Npanels)], ...
        ['contour_SQ_panels' num2str(Npanels)]};
    
    for i=1:6
        print(i,name{i},'-dpng','-r400')
    end
    
end

disp('Done!')

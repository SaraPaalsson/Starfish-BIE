%% Solve Stokes eq. w/ BIE, using specialquad.
close all;
% clear all;
clc

addpath('../mex')

compErrorEst = 0;

% ----------------- Set parameters ----------------------------
res_interf = 'low'; %superlow,low, high
res_domain = 'low'; %superlow, verylow, low, high
interf_param = 'starfish'; %'circle','starfish','ellipse'
typeplot = 'lineplot'; %'filledplot','lineplot'

res = struct(); %result struct
savePlots = 0; %if save plots
savedata = 0;

% ----------------- Setup domain  ------------------------------
[dom] = main_init(res_interf,res_domain,interf_param,typeplot);
res.dom = dom;

% ----------------- Set up problem -----------------------------
% Define boundary condition

% BC caused by point source located at x0 with strength m
x0(1,1) = 2 + 3i;
f0(1,1) = 4*pi + 4*pi*1i;
x0(2,1) = -1.4 - 1.3i;
f0(2,1) = pi - 2*pi*1i;
% % Sum all stokeslets for rhs.
RHS = @(x) comprhs_stokes(x,x0,f0);
res.RHS = RHS;


% ----------------- Calculate density ------------------------------------
% Solve BIE to obtain density mu
mu_stokes = mubie_stokes(dom.N,dom.zDrops,dom.taup(dom.tpar), ...
    dom.taupp(dom.tpar),dom.wDrops,RHS);
res.mu = mu_stokes;

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
evalinterf = 0; %If we compute u on interface
disp('Compute u normal quadrature')
tic
u = compu_stokes(mu_stokes, dom.z, dom.zDrops, ...
    dom.taup(dom.tpar), dom.taupp(dom.tpar), dom.wDrops,evalinterf);
toc

if savedata
    savestr = ['../results/stokes_normq_' res_domain];
    save(savestr,'u')
end
res.u = u;

disp('Compute u special quadrature')
tic
[uspec] = specquad_stokes(u, mu_stokes, dom.Npanels, dom.tau(dom.panels), dom.zDrops, ...
    dom.taup(dom.tpar), dom.wDrops, dom.z, IP1, IP2, W16, W32,u_known);
toc

if savedata
    savestr = ['../../results/stokes_specq_' res_domain];
    save(savestr,'uspec')
end
res.uspec = uspec;

% -------------------------
% Load precomputed resultfiles
%
% load('../../results/normquadu_high')
% load('../../results/specquadu_high')
% load('../../results/errorest_high')
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

disp(['Max err standard quad: e_stand = ' num2str(max(er))])
disp(['Max err special quad: e_spec = ' num2str(max(er_spec))])


error = {reshape(er,size(dom.zplot))  reshape(er_spec,size(dom.zplot))};

if savedata
    savestr = ['../../results/error_' res_domain];
    save(savestr,'error')
end
res.error = error;

%% Compute estimates
compErrorEst = 1
if compErrorEst
    disp('Compute estimates')
    rho1 = 1/(1i*pi)*mu_stokes;
    ntz = -1i*dom.taup(dom.tpar)./abs(dom.taup(dom.tpar));
    rho2 = -1/(2*pi)*mu_stokes.*conj(ntz).^2;
    rho3 = -1/(2*pi)*mu_stokes; %Or should this be mu_stokes*conj(tau-z)?
    
    % Compute error estimate on grid
    errest = error_estL_stokes(dom.z,dom.tau(dom.panels),dom.zDrops,dom.Npanels, ...
        rho1,rho2,rho3,dom.wDrops,dom.taup(dom.tpar));
    
    savestr = ['../../results/errorest_' res_domain];
    save(savestr,'errest')
    res.errest = errest;
end

if strcmp(typeplot,'filledplot')
    levels = -15:3:-3;
    figure(1);
    clf
    contourf(real(dom.zplot),imag(dom.zplot),log10(error{1}),levels,'EdgeColor','none')
    shading flat
    colormap parula
    c = colorbar;
    ylabel(c,'$\log_{10}$ rel. error','FontSize',13,'Interpreter','latex')
    set(gca,'yaxislocation','right');
    set(gca,'YTicklabel',[])
    set(gca,'XTicklabel',[])
    hold on
    plot(real(dom.tau(dom.tpar)), imag(dom.tau(dom.tpar)),'k')
    plot(real(dom.zDrops),imag(dom.zDrops),'.k','MarkerSize',3)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis equal
    % axis([0.5 1.3 0 0.5])
    caxis([-15 0])
    box on
    
    %Add estimate
    errestshape = reshape(errest,size(error{1}));
    contour(real(dom.zplot),imag(dom.zplot),log10(errestshape),levels,'k','linewidth',2)
else
   figure(1); clf
   subplot(121);
   plot(dom.zDrops,'k'); hold on; plot(dom.z,'r*')
   subplot(122);
   loglog(1-dom.reld,(error{1})); hold on; grid on
   loglog(1-dom.reld,errest,'k--')
end
   


%%
%----------------------------------------------------
% Plot
close all

disp('Plot!')
spec_str1 = '10log error normal quad., 35 panels'; spec_str2 = '10log error special quad., 35 panels';
titstr = {spec_str1 spec_str2};

% For contour plot
levels = -15:3:-3;
zoombox = [0.45 1.1 0.1 0.43];
drawbox = @(x,y) plot([x(1) x(2) x(2) x(1) x(1)],[y(1) y(1) y(2) y(2) y(1)],'-k');
tplot = linspace(0,2*pi,1000);
for i=1:2
    figure(i);
    clf
    %     publication_fig
    pcolor(real(dom.zplot),imag(dom.zplot),log10(error{i}))
    shading flat
    % shading interp
    colormap parula
    c = colorbar;
    ylabel(c,'$\log_{10}$ rel. error','FontSize',13,'Interpreter','latex')
    %     caxis([-18 -2])
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
figure(3); clf;
pcolor(real(dom.zplot),imag(dom.zplot),real(RHS(dom.zplot)))
shading flat
hold on
plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
set(gca,'yaxislocation','right');
set(gca,'YTicklabel',[])
set(gca,'XTicklabel',[])
colormap parula
c = colorbar;
ylabel(c,'real$(\mathbf{u})$','FontSize',13,'interpreter','latex')
set(gca,'xtick',[])
set(gca,'ytick',[])
axis equal
% axis([0 1.3 0 1.3])
publication_fig
box on


%%
if 0
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
    % caxis([-15 0])
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

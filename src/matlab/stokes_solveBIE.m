function stokes_solveBIE(varargin)
% Solve Stokes eq. w/ BIE, using specialquad.

% === Set default input values.
compEst = 1;
doSpec = 0;
saveData = 1;
loadPrecom = 0;
for j=1:nargin
    switch j
        case 1
            if ~isempty(varargin{j})
                compEst = varargin{j};
            end
        case 2
            if ~isempty(varargin{j})
                doSpec = varargin{j};
            end
        case 3
            if ~isempty(varargin{j})
                saveData = varargin{j};
            end
        case 4
            if ~isempty(varargin{j})
                loadPrecom = varargin{j};
            end
    end
end


addpath('../mex')

% ----------------- Set parameters ----------------------------
res_interf = 'low'; %superlow,low, high
res_dom = 'low'; %superlow, verylow, low, high
interf_param = 'starfish'; %'circle','starfish','ellipse'
typeplot = 'filledplot';

if ~loadPrecom
    
    res = struct(); %result struct
    
    % ----------------- Setup domain  ------------------------------
    [dom] = main_init(res_interf,res_dom,interf_param,typeplot);
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
    uknown = RHS(dom.z);
    res.uknown = uknown;
    
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
    res.u = u;
    
    % === Special quadrature corrections
    if doSpec
        disp('Compute u special quadrature')
        tic
        [uspec] = specquad_stokes(u, mu_stokes, dom.Npanels, dom.tau(dom.panels), dom.zDrops, ...
            dom.taup(dom.tpar), dom.wDrops, dom.z, IP1, IP2, W16, W32,uknown);
        toc
        
        res.uspec = uspec;
    end
    
    % === Compute errors and plot
    disp('Compute  errors')

    relnorm = norm(uknown,Inf);
    er = abs(uknown-u)/relnorm; er_spec = er;
    Er = reshape(log10(er),size(dom.zplot));
    disp(['Max err standard quad: e_stand = ' num2str(max(er))])
    if doSpec
        er_spec = abs(uknown-uspec)/relnorm;
        Erspec = reshape(log10(er_spec),size(dom.zplot));
        disp(['Max err special quad: e_spec = ' num2str(max(er_spec))])
    end
    error = {reshape(er,size(dom.zplot))  reshape(er_spec,size(dom.zplot))};
    res.error = error;

    % Compute estimates
    if compEst
        disp('Compute estimates')
        tic
        errest = stokes_esterror(dom.z,dom.tau(dom.panels),dom.zDrops,dom.Npanels, ...
            mu_stokes,dom.wDrops,dom.taup(dom.tpar));
        toc
        res.errest = errest;    
    end
    
    if saveData
       disp(['Save data to: results/stokes_D' res_dom '_I' res_interf])
       save(['results/stokes_D' res_dom '_Inps' num2str(dom.Npanels)]); 
    end

    
else
    % === Load precomputed resultfiles
    load(['results/stokes_D' res_dom '_I' res_interf])
    disp(['Using precomputed values from file: results/stokes_D' res_dom '_I' res_interf])
end


% === Plot
figure(1)
clf
pcolor(real(res.dom.zplot), imag(res.dom.zplot), log10(res.error{1}));

num_levels = 7;
levels = linspace(-14, 0, num_levels+1);
colormap('default');
cmap = colormap();
idx = round(linspace(1, size(cmap,1), num_levels));
cmap = cmap( idx, :);
colormap(cmap);

cbar = colorbar();
set(cbar, 'Ticks', levels)
ylabel(cbar,'$\log_{10}$ error','FontSize',20,'interpreter','latex')
shading interp
caxis([min(levels), max(levels)])
hold on
if compEst
    title('Error and estimate - Stokes equation with standard GL quadrature')
    contour(real(res.dom.zplot), imag(res.dom.zplot), log10(abs(res.error{1})), levels(2:end), 'k')
else
    title('Error - Stokes equation with standard GL quadrature')
end
plot(complex([res.dom.zDrops; res.dom.zDrops(1)]), '.-k')
axis equal
axis([-1.4 1.4 -1.4 1.4])
box on
set(gca,'Visible','off')

if doSpec
    figure(2)
    clf
    pcolor(real(res.dom.zplot), imag(res.dom.zplot), log10(res.error{2}));
    title('Error - Stokes equation with special quadrature')
    cbar = colorbar();
    set(cbar, 'Ticks', levels)
    shading interp
    caxis([min(levels), max(levels)])
    hold on
    plot(complex(res.dom.zDrops), '.-k')
    axis equal
    axis([0 1.5 0 1.5])
    box on
    set(gca,'Visible','Off')
end

figure(3); 
clf
uknown = reshape(res.uknown,size(res.dom.zplot));
pcolor(real(res.dom.zplot), imag(res.dom.zplot), abs(uknown));
shading flat; hold on; axis equal
colormap('default');
plot(complex([res.dom.zDrops; res.dom.zDrops(1)]), '.-k')
axis([-1.4 1.4 -1.4 1.4])
h = colorbar;
ylabel(h,'$\mathbf{u}$','interpreter','latex','FontSize',20,'Rotation',0);
set(gca,'Visible','Off')

disp('Done!')
end

% === Extra functions
function u = compu_stokes(mu, z, zDrops, zpDrops, zppDrops, wDrops,evalinterf)
% Compute the velocity field by BIE for stokes equations, for all domain
% points z.
%
% OBS. We assume z not on boundary!
%
% Created 2017-05-28
if evalinterf
    z = zDrops;
    u = zeros(size(z));
    
    M1_jj = imag(0.5*wDrops.*zppDrops./zpDrops);
    M2_jj = 0.5*imag(wDrops.*zppDrops.*conj(zpDrops))./(conj(zpDrops).^2);
    
    for j=1:length(z)
        
        M1 = mu.*wDrops.*imag(zpDrops./(zDrops-z(j)));
        M1(j) = mu(j)*M1_jj(j);
        M2 = conj(mu).*wDrops.*imag(zpDrops.*conj(zDrops-z(j)))./(conj(zDrops-z(j)).^2);
        M2(j) = conj(mu(j))*M2_jj(j);
        
        u(j) = -1i/pi*sum(M1) + 1i/pi*sum(M2);
    end
    
    
else
    u = zeros(size(z));
    
    parfor j=1:length(z)
        %     M1 = sum(mu.*wDrops.*real(zpDrops./(zDrops-z(j))));
        %     M2 = sum(conj(mu).*wDrops.*imag(zpDrops.*conj(zDrops-z(j)))./(conj(zDrops-z(j)).^2));
        %     u(j) = -1/pi*M1 + 1i/pi*M2;
        
        M1 = sum(mu.*wDrops.*imag(zpDrops./(zDrops-z(j))));
        M2 = sum(conj(mu).*wDrops.*imag(zpDrops.*conj(zDrops-z(j)))./(conj(zDrops-z(j)).^2));
        
        u(j) = - 1i/pi*M1 + 1i/pi*M2;
    end
    
end

end

function mu_stokes = mubie_stokes(N,zDrops,zpDrops,zppDrops,wDrops,RHS)
% Calculate density mu from the boundary integral formulation for Stokes
% eq. with Dirichlet boundary condition (RHS)
% Singular integrals calcualted with limit points

% Compute rhs = -h, where h is bc
r = RHS(zDrops);
h1 = -imag(r);
h2 = real(r);
h = h1 + 1i*h2;
rhs = [real(h); imag(h)];

Afunc = @(x) mubie_gmres(x,zDrops,zpDrops,zppDrops,wDrops,[]);

w = gmres(Afunc,rhs,[],1e-13,100);
wr = w(1:N);
wi = w(N+1:end);
mu_stokes = wr + 1i*wi;

end

function RHS = comprhs_stokes(x,x0,f)
% Compute bc for Stokes, using sum of stokeslets

% Number of point forces
Nf = size(f,1);

RHS = 0;
for j=1:Nf
    xhat = x-x0(j);
    r = abs(xhat);
    xh = real(xhat); yh = imag(xhat);
    f1 = real(f(j)); f2 = imag(f(j));
    S11 = -log(r) + (xh.^2)./(r.^2);
    S12 = (xh.*yh)./(r.^2);
    S21 = S12;
    S22 = -log(r) + (yh.^2)./(r.^2);
    RHS1 = S11*f1 + S12*f2;
    RHS2 = S21*f1 + S22*f2;
    RHS = RHS + (RHS1 + 1i*RHS2)/(4*pi);
end

end

% % Extra
% if strcmp(typeplot,'filledplot')
%     levels = -15:3:-3;
%     figure(1);
%     clf
%     err_johan = error{1};
%     err_johan(err_johan < 1e-15) = 1e-15;
%     contourf(real(dom.zplot),imag(dom.zplot),log10(err_johan),levels,'EdgeColor','none')
%     shading flat
%     colormap parula
%     c = colorbar;
%     ylabel(c,'$\log_{10}$ rel. error','FontSize',13,'Interpreter','latex')
%     set(gca,'yaxislocation','right');
%     set(gca,'YTicklabel',[])
%     set(gca,'XTicklabel',[])
%     hold on
%     plot(real(dom.tau(dom.tpar)), imag(dom.tau(dom.tpar)),'k')
%     plot(real(dom.zDrops),imag(dom.zDrops),'.k','MarkerSize',3)
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     axis equal
%     axis([0.5 1.3 0 0.5])
%     caxis([-15 0])
%     box on
%
%     %Add estimate
%     errestshape = reshape(errest,size(error{1}));
%     contour(real(dom.zplot),imag(dom.zplot),log10(errestshape),levels,'k','linewidth',3)
% else
%     figure(1);
%     subplot(121);
%     plot(dom.zDrops,'k'); hold on;
%     cmap = colormap;
%     h1 = plot(dom.z,'*','Color',cmap(54,:))
%     subplot(122);
%     loglog(1-dom.reld,(error{1}),'.','MarkerSize',10,'Color',h1.Color);
%     hold on; grid on
%     loglog(1-dom.reld,errest,'k--','LineWidth',2)
%     ylim([eps(), 10])
%
%
% end
%
%
%
% %%
% %----------------------------------------------------
% % Plot
% close all
%
% disp('Plot!')
% spec_str1 = '10log error normal quad., 35 panels'; spec_str2 = '10log error special quad., 35 panels';
% titstr = {spec_str1 spec_str2};
%
% % For contour plot
% levels = -15:3:-3;
% zoombox = [0.45 1.1 0.1 0.43];
% drawbox = @(x,y) plot([x(1) x(2) x(2) x(1) x(1)],[y(1) y(1) y(2) y(2) y(1)],'-k');
% tplot = linspace(0,2*pi,1000);
% for i=1:2
%     figure(i);
%     clf
%     %     publication_fig
%     pcolor(real(dom.zplot),imag(dom.zplot),log10(error{i}))
%     shading flat
%     % shading interp
%     colormap parula
%     c = colorbar;
%     ylabel(c,'$\log_{10}$ rel. error','FontSize',13,'Interpreter','latex')
%     %     caxis([-18 -2])
%     set(gca,'yaxislocation','right');
%     set(gca,'YTicklabel',[])
%     set(gca,'XTicklabel',[])
%     hold on
%     plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
%     plot(real(dom.zDrops),imag(dom.zDrops),'.k','MarkerSize',3)
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     axis equal
%     axis([0 1.3 0 1.3])
%     caxis([-15 0])
%     publication_fig
%     box on
% end
% figure(3); clf;
% pcolor(real(dom.zplot),imag(dom.zplot),real(RHS(dom.zplot)))
% shading flat
% hold on
% plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
% set(gca,'yaxislocation','right');
% set(gca,'YTicklabel',[])
% set(gca,'XTicklabel',[])
% colormap parula
% c = colorbar;
% ylabel(c,'real$(\mathbf{u})$','FontSize',13,'interpreter','latex')
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% axis equal
% % axis([0 1.3 0 1.3])
% publication_fig
% box on
%
%
% %%
% if 0
%     sfigure(3);
%     clf
%     publication_fig
%     pcolor(real(dom.zplot),imag(dom.zplot),log10(error{2}))
%     shading flat
%     % shading interp
%     colormap parula
%     colorbar
%     caxis([-18 -2])
%     set(gca,'yaxislocation','right');
%     set(gca,'YTicklabel',[])
%     set(gca,'XTicklabel',[])
%     hold on
%     plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
%     set(gca,'xtick',[])
%     set(gca,'ytick',[])
%     axis equal
%     axis([0 1.3 0 1.3])
%     % caxis([-15 0])
%     publication_fig
%     box on
%
%     titstr = {'Level curves, normal quad.' 'Level curves, special quad.'};
%     for i=4:5
%         sfigure(i);
%         clf;
%         publication_fig
%         plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
%         hold on
%         contour(real(dom.zplot),imag(dom.zplot),log10(error{1}),levels,'k')
%         plot(real(dom.zDrops),imag(dom.zDrops),'.k')
%         shading flat
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         axis equal
%         axis(zoombox)
%         axis([0 1.3 0 1.3])
%         caxis([-15 0])
%         axis(zoombox)
%         publication_fig
%         box on
%     end
%
%     for i=6:6
%         sfigure(i);
%         clf
%         publication_fig
%         plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
%         hold on
%         contour(real(dom.zplot),imag(dom.zplot),log10(error{2}),levels,'k')
%         plot(real(dom.zDrops),imag(dom.zDrops),'.k')
%         shading flat
%         set(gca,'xtick',[])
%         set(gca,'ytick',[])
%         axis equal
%         axis(zoombox)
%         axis([0 1.3 0 1.3])
%         caxis([-15 0])
%         axis(zoombox)
%         publication_fig
%         box on
%     end
% end
%
% %----------------------------------------------------
% % Save plots if wanting to
% if savePlots
%     % disp('Save plots')
%     sfigure(1)
%     print -dpng -r400 '../graphics/starfish_error_normq'
%     %crop('../graphics/starfish_error_normq');
%     sfigure(2)
%     print -dpng -r400 '../graphics/starfish_error_interpq'
%     %crop('../graphics/starfish_error_interpq');
%     sfigure(3)
%     print -dpng -r400 '../graphics/starfish_contour_normq'
%     %crop('../graphics/starfish_contour_normq');
%     sfigure(4)
%     print -dpng -r400 '../graphics/starfish_contour_interpq'
%     %crop('../graphics/starfish_contour_interpq');
%
%     name = {['filled_error_panels' num2str(Npanels)], ...
%         ['fillederror_SQbox_panels' num2str(Npanels)], ...
%         ['fillederror_SQ_panels' num2str(Npanels)], ...
%         ['contour_panels' num2str(Npanels)], ...
%         ['contour__LC_panels' num2str(Npanels)], ...
%         ['contour_SQ_panels' num2str(Npanels)]};
%
%     for i=1:6
%         print(i,name{i},'-dpng','-r400')
%     end
%
% end

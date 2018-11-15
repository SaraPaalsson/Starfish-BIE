function stokes_solveBIE(varargin)
% Solve Stokes eq. w/ BIE, using specialquad.
varnms={'compEst','doSpec','saveData','loadPrecom','preComFile','fileOut'};
defvals={'0','1','0','0','''''','''tmp'''};

%assign default values:
for i=1:length(varnms)
    eval([varnms{i} '=' defvals{i} ';']);
end
%assign specific values
for i= 1:nargin
    if ischar(varargin{i})
    varargin{i}=['''' varargin{i} ''''];    
    else
    varargin{i}=num2str(varargin{i});     
    end
    
    if strcmpi(varnms{i},'~')==false
            eval([varnms{i} '=' varargin{i} ';']);
    end
end
   

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
    x0(1,1) = 1.1 + 1.3i;
    f0(1,1) = 4*pi + 4*pi*1i;
    x0(2,1) = -1.4 - 1.3i;
    f0(2,1) = pi/ - 2*pi*1i;
    x0(3,1) = 1.3-0.75i;
    f0(3,1) = -0.5*pi+3.5*pi*1i;
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
%        save(['results/stokes_D' res_dom '_Inps' num2str(dom.Npanels)]);
        save(['results/' fileOut '_D' res_dom '_Inps' num2str(dom.Npanels)]);
    end

    
else
    % === Load precomputed resultfiles
    load(['results/' preComFile])
%     load(['results/stokes_D' res_dom '_I' res_interf])
    disp(['Using precomputed values from file: results/stokes_D' res_dom '_I' res_interf])

    disp(['Compute estimates: ' num2str(compEst)])
    disp(['Compute special quadrature: ' num2str(doSpec)])
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
set(cbar, 'Ticks', levels,'FontSize',15)
ylabel(cbar,'$\log_{10}$ error','FontSize',25,'interpreter','latex')
shading interp
caxis([min(levels), max(levels)])
hold on
if compEst
%     title('Error and estimate - Stokes equation with standard GL quadrature')
    contour(real(res.dom.zplot), imag(res.dom.zplot), log10(abs(res.error{1})), levels(2:end), 'k')
else
%     title('Error - Stokes equation with standard GL quadrature')
end
plot(complex([res.dom.zDrops; res.dom.zDrops(1)]), '-k','LineWidth',2)
axis equal
box on

axis([0 1.35 0 1.35])

set(gca,'Visible','off')

% % % Make box
xmin = 0.35; xmax = 1.35;
ymin = 0; ymax = 0.5;
% Plot box
plot([xmin xmax],[ymin ymin],'k-','LineWidth',2)
plot([xmin xmax],[ymax ymax],'k-','LineWidth',2)
plot([xmin xmin],[ymin ymax],'k-','LineWidth',2)
plot([xmax xmax],[ymin ymax],'k-','LineWidth',2)
% Cut box
axis([xmin xmax ymin ymax])

% set(cbar,'Visible','Off')

if doSpec
    figure(2)
    clf
    pcolor(real(res.dom.zplot), imag(res.dom.zplot), log10(res.error{2}));
    title('Error - Stokes equation with special quadrature')
    cbar = colorbar();
    set(cbar, 'Ticks', levels,'FontSize',15)
    shading interp
    ylabel(cbar,'$\log_{10}$ error','FontSize',25,'interpreter','latex')

    num_levels = 7;
    levels = linspace(-14, 0, num_levels+1);
%     colormap('default');
    cmap = colormap();
    idx = round(linspace(1, size(cmap,1), num_levels));
    cmap = cmap( idx, :);
    colormap(cmap);

    
    caxis([min(levels), max(levels)])
    hold on
    plot(complex(res.dom.zDrops), '-k','LineWidth',2,'MarkerSize',10)
    axis equal
%     axis([0 1.5 0 1.5])
    axis([-1.2 1.4 -1.4 1.4])
    box on
    set(gca,'Visible','Off')
    
% % % % Make box
% xmin = 0.35; xmax = 1.35;
% ymin = 0; ymax = 0.5;
% % % Quadrant
xmin = 0; xmax = 1.35;
ymin = 0; ymax = 1.25;

%     % % Plot box
%     plot([xmin xmax],[ymin ymin],'k-','LineWidth',2)
%     plot([xmin xmax],[ymax ymax],'k-','LineWidth',2)
%     plot([xmin xmin],[ymin ymax],'k-','LineWidth',2)
%     plot([xmax xmax],[ymin ymax],'k-','LineWidth',2)
    % Cut box
%     axis([xmin xmax ymin ymax])


    set(cbar,'Visible','Off')

end

figure(3); 
clf
uknown = reshape(res.uknown,size(res.dom.zplot));
pcolor(real(res.dom.zplot), imag(res.dom.zplot), abs(uknown));
shading flat; hold on; axis equal
colormap('default');
plot(complex([res.dom.zDrops; res.dom.zDrops(1)]), '-k','LineWidth',2)
% Plot point sources
plot(x0(1),'k.','MarkerSize',35)
plot(x0(2),'k.','MarkerSize',35)
plot(x0(3),'k.','MarkerSize',35)
axis([-2 2 -2 2])
box on
xlabel('$\Re(z)$','FontSize',25)
ylabel('$\Im(z)$','FontSize',25)
set(gca,'FontSize',20)
h = colorbar;
set(h,'FontSize',15)
ylabel(h,'$|u|$','interpreter','latex','FontSize',25,'Rotation',0);
% set(gca,'Visible','Off')


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

Afunc = @(x) mubie_gmres(x,zDrops,zpDrops,zppDrops,wDrops);

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

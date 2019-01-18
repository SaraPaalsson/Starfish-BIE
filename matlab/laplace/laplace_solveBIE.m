function laplace_solveBIE(varargin)
% Solve BIE for starfish by computing complex density. Correct for nearly
% singular integrals using special quadrature. Compute error estimates.


% === Set default input values.
compEst = 0;
doSpec = 1;
saveData = 0;
loadPrecom = 0;
fileN = '';
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
        case 5
            if ~isempty(varargin{j})
                fileN = varargin{j};
            end
    end
end
addpath('../mex')


% === Set parameters
res_interf = 'low'; %superlow,low, high
res_domain = 'low'; %superlow, verylow, low, high
interf_param = 'starfish';
typeplot = 'filledplot'; %'filledplot','lineplot'


if strcmp(fileN,'')
    fileN = ['_D' res_domain '_I' res_interf];
end

% === Setup domain
res = struct(); %result struct
dom = main_init(res_interf,res_domain,interf_param,typeplot);
res.dom = dom;


% === Set up problem
% Define boundary condition, Laplace's eq. --> Provides exact solution!
zsrc1 = 3+3i;
zsrc2 = -2.5-2.5i;
RHS = @(x) real( 1 ./ (x-zsrc1) + 1 ./ (x-zsrc2));
res.RHS = RHS;


% === Calculate density
% Solve BIE to obtain density mu
mu_lapl = mubie_lapl(dom.N,dom.zDrops,dom.taup(dom.tpar), ...
    dom.taupp(dom.tpar),dom.wDrops,RHS);
res.mu = mu_lapl;


% === Calculate known solution over the domain
u_known = RHS(dom.z);
res.uknown = u_known;

if ~loadPrecom
    % === Compute u with normal and special quadrature
    load 'IP1632.mat'
    load 'glW.mat' %read in GL 16 and 32 weights
    
    
    % === Compute u over the domain
    % Use 16-GL when possible
    % Use special quadrature for points too close to the boundary
    disp('Compute u normal quadrature')
    tic
    u = laplace_compu(mu_lapl, dom.z, dom.zDrops, dom.taup(dom.tpar), dom.wDrops);
    toc
    
    if saveData
        savestr = ['results/normquadu' fileN];
        %         savestr = ['results/normquadu_D' res_domain '_I' res_interf];
        save(savestr,'u')
    end
    res.u = u;
    
    if doSpec
        disp('Compute u special quadrature')
        tic
        [uspec] = laplace_specquad(u, mu_lapl, dom.Npanels, dom.tau(dom.panels), dom.zDrops, ...
            dom.taup(dom.tpar), dom.wDrops, dom.z, IP1, IP2, W16, W32);
        toc
        % disp('Compute u special quadrature MEX')
        % tic
        % [uspec,~] = mex_saraspecquad(u, mu_lapl, dom.tau(dom.panels), dom.zDrops, dom.taup(dom.tpar), dom.wDrops, dom.z);
        % toc
        
        if saveData
            savestr = ['results/specquadu' fileN];
            %             savestr = ['results/specquadu_D' res_domain '_I' res_interf];
            save(savestr,'uspec')
        end
        res.uspec = uspec;
    end
    
    % === Compute error estimates
    if compEst
        disp('Compute estimates')
        tic
        errest = laplace_esterror(dom.z,dom.tau(dom.panels),dom.zDrops,dom.Npanels, ...
            mu_lapl,dom.wDrops,dom.taup(dom.tpar));
        toc
        if saveData
            savestr = ['results/errorest' fileN];
            %             savestr = ['results/errorest_D' res_domain '_I' res_interf];
            save(savestr,'errest')
        end
        res.errest = errest;
        est = reshape(errest,size(dom.zplot));
    end
else
    % === Load precomputed resultfiles
    load(['results/normquadu' fileN])
    load(['results/specquadu' fileN])
    load(['results/errorest' fileN])
    
%     load(['results/normquadu_D' res_domain '_I' res_interf])
%     load(['results/specquadu_D' res_domain '_I' res_interf])
%     load(['results/errorest_D' res_domain '_I' res_interf])
end


% === Compute errors and plot
disp('Compute  errors')

relnorm = norm(u_known,Inf);
er = abs(u_known-u)/relnorm; er_spec = er;
Er = reshape(log10(er),size(dom.zplot));
if doSpec
    er_spec = abs(u_known-uspec)/relnorm;
    Erspec = reshape(log10(er_spec),size(dom.zplot));
end

error = {reshape(er,size(dom.zplot))  reshape(er_spec,size(dom.zplot))};

if saveData
    savestr = ['../../results/error_' res_domain];
    save(savestr,'error')
end
res.error = error;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% === Plot
figure(1)
clf
pcolor(real(dom.zplot), imag(dom.zplot), log10(error{1}+eps));
num_levels = 7;
levels = linspace(-14, 0, num_levels+1);
colormap('default');
cmap = colormap();
idx = round(linspace(1, size(cmap,1), num_levels));
cmap = cmap( idx, :);
colormap(cmap);

cbar = colorbar();
set(cbar, 'Ticks', levels)
hh=ylabel(cbar,'$\log_{10}$-error');
set(hh,'interpreter','latex','FontSize',20)
shading interp
caxis([min(levels), max(levels)])
hold on
if compEst
    contour(real(dom.zplot), imag(dom.zplot), log10(abs(error{1})), levels(2:end), 'k')
end
plot(complex(dom.zDrops), 'k-')
axis equal
% axis([0 1.5 0 1.5])
ylim([-1.3 1.3])
% set(gca,'Visible','Off')
box on

% Draw box
xmax = 1.4;
xmin = 0.3;
ymax = 0.5;
ymin = 0;
plot([xmin xmax; xmin xmax; xmin xmin; xmax xmax],[ymin ymin; ymax ymax; ymin ymax; ymin ymax],'k-')


% === Plot
figure(11)
clf
pcolor(real(dom.zplot), imag(dom.zplot), log10(error{1}+eps));
num_levels = 7;
levels = linspace(-14, 0, num_levels+1);
colormap('default');
cmap = colormap();
idx = round(linspace(1, size(cmap,1), num_levels));
cmap = cmap( idx, :);
colormap(cmap);

cbar = colorbar();
set(cbar, 'Ticks', levels)
hh=ylabel(cbar,'$\log_{10}$-error');
set(hh,'interpreter','latex','FontSize',20)
shading interp
caxis([min(levels), max(levels)])
hold on
if compEst
    contour(real(dom.zplot), imag(dom.zplot), log10(abs(error{1})), levels(2:end), 'k')
end
plot(complex(dom.zDrops), 'k-')
axis equal
axis([xmin xmax ymin ymax])
box on

% Draw box
xmax = 1.4;
xmin = 0.3;
ymax = 0.5;
ymin = 0;
plot([xmin xmax; xmin xmax; xmin xmin; xmax xmax],[ymin ymin; ymax ymax; ymin ymax; ymin ymax],'k-')


if doSpec
    figure(2)
    clf
    pcolor(real(dom.zplot), imag(dom.zplot), log10(error{2}+eps));
    %     title('Error - Laplace equation with special quadrature')
    cbar = colorbar();
    set(cbar, 'Ticks', levels)
    hh=ylabel(cbar,'$\log_{10}$-error');
    set(hh,'interpreter','latex','FontSize',20)
    shading interp
    caxis([min(levels), max(levels)])
    hold on
    plot(complex(dom.zDrops), 'k-')
    axis equal
    axis([xmin xmax ymin ymax])
    box on
end



disp('Done!')
end

% ------------------------------------------------------------------------
% === ADDITIONAL FUNCTIONS
% ------------------------------------------------------------------------
function laplace_density = mubie_lapl(N,zDrops,zpDrops,zppDrops,wDrops,RHS)
% Compute density mu (complex) from the boundary integral formulation for Laplace
% eq. with Dirichlet boundary condition (RHS)
% Singular integrals calcualted with limit points
A = eye(N,N);
for i=1:N
    A(i,:) = A(i,:) + (1/pi*wDrops.*imag(zpDrops./(zDrops-zDrops(i))))';
end
d = 1 + 1/pi*wDrops.*imag(zppDrops./(2*zpDrops));
b = 2*RHS(zDrops);
A(logical(eye(size(A)))) = 0;
A = A + diag(d,0);
laplace_density = gmres(A,b,[],1e-13,N);
end

function u = laplace_compu(mu_lapl, z, zDrops, zpDrops, wDrops)
% Compute the solution to Laplace's equation given the complex density mu.
u = zeros(size(z));

parfor j=1:length(z)
    u(j) = sum(mu_lapl.*wDrops.*imag(zpDrops./(zDrops-z(j))));
end
u = 1/(2*pi)*u;

end


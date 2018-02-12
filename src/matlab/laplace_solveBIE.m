function laplace_solveBIE(varargin)
% Solve BIE for starfish by computing complex density. Correct for nearly
% singular integrals using special quadrature. Compute error estimates.


% === Set default input values.
compEst = 1;
doSpec = 0;
saveData = 0;
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


if ~loadPrecom
    % === Set parameters
    res_interf = 'low'; %superlow,low, high
    res_domain = 'verylow'; %superlow, verylow, low, high
    interf_param = 'starfish';
    typeplot = 'filledplot'; %'filledplot','lineplot'
    
    
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
    
    % === Compute u with normal and special quadrature
    load 'IP1632.mat'
    load 'glW.mat' %read in GL 16 and 32 weights
    
    
    % === Compute u over the domain
    % Use 16-GL when possible
    % Use special quadrature for points too close to the boundary
    disp('Compute u normal quadrature')
    tic
    u = laplace_compu(dom.N, mu_lapl, dom.z, dom.zDrops, dom.taup(dom.tpar), dom.wDrops);
    toc
    
    if saveData
        savestr = ['../results/normquadu_' res_domain];
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
            savestr = ['../../results/specquadu_' res_domain];
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
            savestr = ['../../results/errorest_' res_domain];
            save(savestr,'errest')
        end
        res.errest = errest;
        est = reshape(errest,size(dom.zplot));
    end
else
    % === Load precomputed resultfiles
    
    load('../../results/normquadu_high')
    load('../../results/specquadu_high')
    load('../../results/errorest_high')
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


% === Plot
figure(1)
clf
pcolor(real(dom.zplot), imag(dom.zplot), log10(error{1}));

num_levels = 7;
levels = linspace(-14, 0, num_levels+1);
colormap('default');
cmap = colormap();
idx = round(linspace(1, size(cmap,1), num_levels));
cmap = cmap( idx, :);
colormap(cmap);

cbar = colorbar();
set(cbar, 'Ticks', levels)
shading interp
caxis([min(levels), max(levels)])
hold on
if compEst
    title('Error and estimate - Laplace equation with standard GL quadrature')
    contour(real(dom.zplot), imag(dom.zplot), log10(abs(error{1})), levels(2:end), 'k')
else
    title('Error - Laplace equation with standard GL quadrature')
end
plot(complex(dom.zDrops), '.-r')
axis equal
axis([0 1.5 0 1.5])
box on

if doSpec
    figure(2)
    clf
    pcolor(real(dom.zplot), imag(dom.zplot), log10(error{2}));
    title('Error - Laplace equation with special quadrature')
    cbar = colorbar();
    set(cbar, 'Ticks', levels)
    shading interp
    caxis([min(levels), max(levels)])
    hold on
    plot(complex(dom.zDrops), '.-r')
    axis equal
    axis([0 1.5 0 1.5])
    box on
end
    


disp('Done!')
end

% === ADDITIONAL FUNCTIONS
function laplace_density = mubie_lapl(N,zDrops,zpDrops,zppDrops,wDrops,RHS)
% Calculate density mu from the boundary integral formulation for Laplace
% eq. with Dirichlet boundary condition (RHS)
% Singular integrals calcualted with limit points
A = eye(N,N);
for i=1:N
    A(i,:) = A(i,:) + (1/pi*wDrops.*imag(zpDrops./(zDrops-zDrops(i))))';
end
d = 1 + 1/pi*wDrops.*imag(zppDrops./(2*zpDrops));
b = 2*RHS(zDrops);%imag(tauk.^2./(tauk-zp));
A(logical(eye(size(A)))) = 0;
A = A + diag(d,0);
laplace_density = gmres(A,b,[],1e-13,N);
% mu_lapl = A\b;
end

function u = laplace_compu(N, mu_lapl, z, zDrops, zpDrops, wDrops)

u = zeros(size(z));

parfor j=1:length(z)
    u(j) = sum(mu_lapl.*wDrops.*imag(zpDrops./(zDrops-z(j))));
end
u = 1/(2*pi)*u;

end


% % % % Addition plot things
% % % 
% % % %----------------------------------------------------
% % % % Plot
% 
% figure(1);
% clf
% contour(real(dom.zplot),imag(dom.zplot),log10(error{1}),levels,'k')
% % pcolor(real(dom.zplot),imag(dom.zplot),log10(error{1}))
% shading flat
% % shading interp
% colormap parula
% colorbar
% set(gca,'yaxislocation','right');
% set(gca,'YTicklabel',[])
% set(gca,'XTicklabel',[])
% hold on
% % plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
% plot(real(dom.zDrops),imag(dom.zDrops),'.-k','MarkerSize',10)
% set(gca,'xtick',[])
% set(gca,'ytick',[])
% % axis equal
% caxis([-15 0])
% % publication_fig
% box on
% axis square
% axis([0 2.1 0 2.1])
% 
% 
% 
% % % disp('Plot!')
% % % switch typeplot
% % %     case 'filledplot'
% % %         spec_str1 = '10log error normal quad., 35 panels'; spec_str2 = '10log error special quad., 35 panels';
% % %         titstr = {spec_str1 spec_str2};
% % %         
% % %         % For contour plot
% % %         levels = -15:3:-3;
% % %         zoombox = [0.45 1.1 0.1 0.43];
% % %         drawbox = @(x,y) plot([x(1) x(2) x(2) x(1) x(1)],[y(1) y(1) y(2) y(2) y(1)],'-k');
% % %         tplot = linspace(0,2*pi,1000);
% % %         
% % %         for i=1:2
% % %             sfigure(i);
% % %             clf
% % %             publication_fig
% % %             pcolor(real(dom.zplot),imag(dom.zplot),log10(error{i}))
% % %             shading flat
% % %             % shading interp
% % %             colormap parula
% % %             colorbar
% % %             caxis([-18 -2])
% % %             set(gca,'yaxislocation','right');
% % %             set(gca,'YTicklabel',[])
% % %             set(gca,'XTicklabel',[])
% % %             hold on
% % %             plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
% % %             plot(real(dom.zDrops),imag(dom.zDrops),'.k','MarkerSize',3)
% % %             set(gca,'xtick',[])
% % %             set(gca,'ytick',[])
% % %             axis equal
% % %             %             axis([0 1.3 0 1.3])
% % %             caxis([-15 0])
% % %             publication_fig
% % %             box on
% % %         end
% % %         
% % %         if 1
% % %             sfigure(3);
% % %             clf
% % %             publication_fig
% % %             pcolor(real(dom.zplot),imag(dom.zplot),log10(error{2}))
% % %             shading flat
% % %             % shading interp
% % %             colormap parula
% % %             colorbar
% % %             caxis([-18 -2])
% % %             set(gca,'yaxislocation','right');
% % %             set(gca,'YTicklabel',[])
% % %             set(gca,'XTicklabel',[])
% % %             hold on
% % %             plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
% % %             set(gca,'xtick',[])
% % %             set(gca,'ytick',[])
% % %             axis equal
% % %             axis([0 1.3 0 1.3])
% % %             caxis([-15 0])
% % %             publication_fig
% % %             box on
% % %             
% % %             titstr = {'Level curves, normal quad.' 'Level curves, special quad.'};
% % %             for i=4:5
% % %                 sfigure(i);
% % %                 clf;
% % %                 publication_fig
% % %                 plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
% % %                 hold on
% % %                 contour(real(dom.zplot),imag(dom.zplot),log10(error{1}),levels,'k')
% % %                 plot(real(dom.zDrops),imag(dom.zDrops),'.k')
% % %                 shading flat
% % %                 set(gca,'xtick',[])
% % %                 set(gca,'ytick',[])
% % %                 axis equal
% % %                 axis(zoombox)
% % %                 axis([0 1.3 0 1.3])
% % %                 caxis([-15 0])
% % %                 axis(zoombox)
% % %                 publication_fig
% % %                 box on
% % %             end
% % %             
% % %             for i=6:6
% % %                 sfigure(i);
% % %                 clf
% % %                 publication_fig
% % %                 plot(real(dom.tau(tplot)), imag(dom.tau(tplot)),'k')
% % %                 hold on
% % %                 contour(real(dom.zplot),imag(dom.zplot),log10(error{2}),levels,'k')
% % %                 plot(real(dom.zDrops),imag(dom.zDrops),'.k')
% % %                 shading flat
% % %                 set(gca,'xtick',[])
% % %                 set(gca,'ytick',[])
% % %                 axis equal
% % %                 axis(zoombox)
% % %                 axis([0 1.3 0 1.3])
% % %                 caxis([-15 0])
% % %                 axis(zoombox)
% % %                 publication_fig
% % %                 box on
% % %             end
% % %         end
% % %         
% % %         %----------------------------------------------------
% % %         % Save plots if wanting to
% % %         if savePlots
% % %             % disp('Save plots')
% % %             sfigure(1)
% % %             print -dpng -r400 '../graphics/starfish_error_normq'
% % %             %crop('../graphics/starfish_error_normq');
% % %             sfigure(2)
% % %             print -dpng -r400 '../graphics/starfish_error_interpq'
% % %             %crop('../graphics/starfish_error_interpq');
% % %             sfigure(3)
% % %             print -dpng -r400 '../graphics/starfish_contour_normq'
% % %             %crop('../graphics/starfish_contour_normq');
% % %             sfigure(4)
% % %             print -dpng -r400 '../graphics/starfish_contour_interpq'
% % %             %crop('../graphics/starfish_contour_interpq');
% % %             
% % %             name = {['filled_error_panels' num2str(Npanels)], ...
% % %                 ['fillederror_SQbox_panels' num2str(Npanels)], ...
% % %                 ['fillederror_SQ_panels' num2str(Npanels)], ...
% % %                 ['contour_panels' num2str(Npanels)], ...
% % %                 ['contour__LC_panels' num2str(Npanels)], ...
% % %                 ['contour_SQ_panels' num2str(Npanels)]};
% % %             
% % %             for i=1:6
% % %                 print(i,name{i},'-dpng','-r400')
% % %             end
% % %             
% % %         end
% % %         
% % %     case 'lineplot'
% % %         figure();
% % %         semilogy(flipud(dom.reld),error{1},'--','MarkerSize',15,'DisplayName','Standard quad., 35 panels'); hold on;
% % %         semilogy(flipud(dom.reld),error{2},'.-','DisplayName','Special quad., 35 panels')
% % %         legend('toggle')
% % %         grid on
% % %         xlabel('Relative distance to interface, $r$')
% % %         ylabel('Relative error')
% % % end

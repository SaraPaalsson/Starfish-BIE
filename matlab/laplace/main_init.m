function dom  = main_init(res_interf,res_domain,interf_param,typeplot)

% Set number of panels on interface
switch res_interf
    case 'low'
        Npanels = 50;
    case 'high'
        Npanels = 100;
    case 'superlow'
        Npanels = 25;
    case 'superhigh'
        Npanels = 160;
end
N = Npanels*16;

% Create parametrization
switch interf_param
    case 'starfish'
%         s = pi/5;
        s = 0;
%         a = 0.225;
        a = 0.3;
        b = 5;
        tau = @(t) (1+a*cos(b*(t+s))).*exp(1i*(t+s)); %starfish parametrization
        taup = @(t) (-a*b*sin(b*(t+s))+1i*(1+a*cos(b*(t+s)))).*exp(1i*(t+s));
        taupp = @(t) exp(1i*(t+s)).*(-1+(-a*b^2-a)*cos(b*(t+s))-2*a*b*1i*sin(b*(t+s)));
    case 'circle'
        r = 2;
        tau = @(t) r*(cos(t) + 1i*sin(t));
        taup = @(t) r*(-sin(t) + 1i*cos(t));
        taupp = @(t) r*(-cos(t) - 1i*sin(t));
    case 'ellipse'
        a = sqrt(2); b = 0.4;
        tau = @(t) exp(1i*0.4)*(a*cos(t) + 1i*b*sin(t));
        taup = @(t) exp(1i*0.4)*(-a*sin(t) + 1i*b*cos(t));
        taupp = @(t) exp(1i*0.4)*(-a*cos(t) - 1i*b*sin(t));
end

% Obtan GL-16 nodes and weights
[zDrops, wDrops,tpar,panels] = gl16(Npanels,tau);

% Fill the domain with computational points, discretized in r and t
%  --> Grid coarse far from and fine close to boundary
%  --> Compute only on the first quadrant (for speed)

switch typeplot
    case 'filledplot'
        switch res_domain
            case 'superlow'
                nbrR = 40;
                nbrT = 40;
            case 'verylow'
                nbrR = 100;
                nbrT = 100;
            case 'low'
                nbrR = 300;
                nbrT = 300;
            case 'medium'
                nbrR = 500;
                nbrT = 500;
            case 'high'
                nbrR = 2000;
                nbrT = 2500;
        end
        
        % Where to go from coarse to fine grid
        R1 = 0.4;
        r1 = linspace(0,R1,20); %5
        
        r2 = linspace(R1, 0.999,nbrR); r2  = r2(2:end);
        r = [r1 r2]';
                
        t = linspace(0,2*pi,nbrT)';
        t = t(t <= pi/2);

%         t = linspace(0,2*pi,nbrT+1)'; t = t(1:end-1);
        [Rplot,Tplot] = meshgrid(r,t);
        z = Rplot(:).*tau(Tplot(:));
        zplot = Rplot.*tau(Tplot);
        
        dom = struct('z',z,'zplot',zplot,'zDrops',zDrops,'wDrops',wDrops,...
    'tpar',tpar,'panels',panels,'N',N,'Npanels',Npanels,'tau',tau,'taup',taup,'taupp',taupp);

        
    case 'lineplot'
        switch res_domain
            case 'superlow'
                nbrP = 10;
            case 'verylow'
                nbrP = 100;
            case 'low'
                nbrP = 400;
            case 'high'
                nbrP = 2000;
        end
        reld = linspace(0.5,0.99,nbrP)';
        
        
        z = zDrops(end-50)*reld;
        
        
        zplot = z;
                
        dom = struct('z',z,'zplot',zplot,'zDrops',zDrops,'wDrops',wDrops,...
    'tpar',tpar,'panels',panels,'N',N,'Npanels',Npanels,'tau',tau,'taup',taup,'taupp',taupp, ...
    'reld',reld);

        
end


end

function [zDrops, wDrops, tpar,panels] = gl16(Npanels,tau)
% Divide domain parametrized by tau into Npanels panels and distribute 16
% GL nodes and weights on each
%   zDrops, GL nodes in complex plane
%   spar, GL nodes in parametrization [0,2pi]
%   wDrops, GL weights

panels = linspace(0,2*pi,Npanels+1)'; %divide parametrization into panels

zDrops = zeros(Npanels*16,1); %complex interface coord.
tpar = zeros(Npanels*16,1); %parametrization vector
wDrops = zeros(Npanels*16,1); %GL. weights

for i=1:Npanels
    [nodes,weights] = gaussleg(16,[panels(i) panels(i+1)]);
    tpar((i-1)*16+1:16*i) = nodes;
    zDrops((i-1)*16+1:16*i) = tau(nodes);
    wDrops((i-1)*16+1:16*i) = weights;
end

end

function [nodes,weights] = gaussleg(n,interv)
% Creates G-L nodes and weights for polynomials of degree n, on
% the interval interv (as done in Trefethen)
n_vec = 1:(n-1);
beta = 0.5*(1-(2*n_vec).^(-2)).^(-1/2);
Tn = diag(beta,-1) + diag(beta,1);
[V,D] = eig(Tn);
nodes = D(logical(eye(size(D))));
weights = (2*V(1,:).^2)'; %here we could use our saved W16 weights and
%remap them instead... but we won't save much time doing so

nodes = (interv(1)*(1-nodes)+interv(2)*(1+nodes))/2;
weights =(interv(2)-interv(1))/2*weights;
end
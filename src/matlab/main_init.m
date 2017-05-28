function dom  = main_init(res_interf,res_domain,interf_param)

% Set number of panels on interface
switch res_interf
    case 'low'
        Npanels = 35;
    case 'high'
        Npanels = 70;
    case 'superlow'
        Npanels = 20;
end
N = Npanels*16;

% Create parametrization
switch interf_param
    case 'starfish'
        s = 0;
        tau = @(t) (1+0.3*cos(5*(t+s))).*exp(1i*(t+s)); %starfish parametrization
        taup = @(t) (-1.5*sin(5*(t+s))+1i*(1+0.3*cos(5*(t+s)))).*exp(1i*(t+s));
        taupp = @(t) exp(1i*(t+s)).*(-1-7.8*cos(5*(t+s))-(3i)*sin(5*(t+s)));
    case 'circle'
        tau = @(t) cos(t) + 1i*sin(t);
        taup = @(t) -sin(t) + 1i*cos(t);
        taupp = @(t) -cos(t) - 1i*sin(t);
end

% Obtan GL-16 nodes and weights
[zDrops, wDrops,tpar,panels] = gl16(Npanels,tau);

% Fill the domain with computational points, discretized in r and t
%  --> Grid coarse far from and fine close to boundary
%  --> Compute only on the first quadrant (for speed)

switch res_domain
    case 'superlow'
        nbrR = 10;
        nbrT = 10;
    case 'verylow'
        nbrR = 100;
        nbrT = 100;
    case 'low'
        nbrR = 300;
        nbrT = 300;
    case 'high'
        nbrR = 2000;
        nbrT = 2500;
end

% Where to go from coarse to fine grid
R1 = 0.4; 

r = [linspace(0,R1,10) linspace(R1, 0.999,nbrR)]';
t = linspace(0,pi/2,nbrT)';
[Rplot,Tplot] = meshgrid(r,t);
z = Rplot(:).*tau(Tplot(:));
zplot = Rplot.*tau(Tplot);

dom = struct('z',z,'zplot',zplot,'zDrops',zDrops,'wDrops',wDrops,...
    'tpar',tpar,'panels',panels,'N',N,'Npanels',Npanels,'tau',tau,'taup',taup,'taupp',taupp);

end
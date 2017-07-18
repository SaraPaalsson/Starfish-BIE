function dom  = main_init(res_interf,res_domain,interf_param,typeplot)

% Set number of panels on interface
switch res_interf
    case 'low'
        Npanels = 40;
    case 'high'
        Npanels = 80;
    case 'superlow'
        Npanels = 20;
    case 'superhigh'
        Npanels = 160;
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
        % t = linspace(0,pi/2,nbrT)';
        t = linspace(0,2*pi,nbrT)';
        % t = 0;
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
                nbrP = 300;
            case 'high'
                nbrP = 2000;
        end
        reld = linspace(0,0.999,nbrP)';
        z = zDrops(1)*reld;
        zplot = z;
        
        dom = struct('z',z,'zplot',zplot,'zDrops',zDrops,'wDrops',wDrops,...
    'tpar',tpar,'panels',panels,'N',N,'Npanels',Npanels,'tau',tau,'taup',taup,'taupp',taupp, ...
    'reld',reld);

        
end


end
function dom  = main_init(res_interf,res_domain,interf_param,typeplot)

% Set number of panels on interface
switch res_interf
    case 'low'
        Npanels = 24;
    case 'high'
        Npanels = 50;
    case 'superlow'
        Npanels = 12;
    case 'superhigh'
        Npanels = 100;
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
    case 'corner'
        theta =  7 * pi / 4;
        tau = @(t) corner_z(t, theta);
        taup = @(t) corner_zp(t, theta);
        taupp = @(t) corner_zpp(t, theta);
    case 'square'        
        tau = @(t) square_z(t);
        taup = @(t) square_zp(t);
        taupp = @(t) zeros(size(t));
end

% Obtan GL-16 nodes and weights
if strcmp(interf_param, 'corner')
    panels = linspace(-0.5, 0.5,Npanels+1)'; %divide parametrization into panels
else
    panels = linspace(0,2*pi,Npanels+1)'; %divide parametrization into panels
end

[zDrops, wDrops,tpar,panels] = gl16(Npanels,tau, panels);

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
        
        if strcmp(interf_param, 'square') || strcmp(interf_param, 'corner')
            x = linspace(min(real(zDrops)), max(real(zDrops)), nbrT);
            y = linspace(min(imag(zDrops)), max(imag(zDrops)), nbrT);
            [X, Y] = meshgrid(x,y);
            
            zplot = X + 1i * Y;
            z = zplot(:);
        else
            t = linspace(0,2*pi,nbrT)';
            t = t(t <= pi/2);
            
            [Rplot,Tplot] = meshgrid(r,t);
            z = Rplot(:).*tau(Tplot(:));
            zplot = Rplot.*tau(Tplot);
        end
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

function [zDrops, wDrops, tpar,panels] = gl16(Npanels,tau, panels)
% Divide domain parametrized by tau into Npanels panels and distribute 16
% GL nodes and weights on each
%   zDrops, GL nodes in complex plane
%   spar, GL nodes in parametrization [0,2pi]
%   wDrops, GL weights

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = corner_z(s, theta)

x = zeros(size(s));
y = zeros(size(s));

for i = 1:length(s)
    if s(i) < 0
        x(i) = -sin(pi*s(i)).*cos((s(i) + 0.5)*theta);
        y(i) = -sin(pi*s(i)).*sin((s(i) + 0.5)*theta);
    else
        x(i) = sin(pi*s(i)).*cos((s(i) - 0.5)*theta);
        y(i) = sin(pi*s(i)).*sin((s(i) - 0.5)*theta);
    end
end

z = x + 1i*y;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zp = corner_zp(s, theta)

xp = zeros(size(s));
yp = zeros(size(s));

for i = 1:length(s)
    
    if s(i) < 0
        xp(i) = -pi*cos(pi*s(i)).*cos((s(i) + 0.5)*theta) + ...
                    theta * sin(pi*s(i)).*sin((s(i) + 0.5)*theta);
        yp(i) = -pi*cos(pi*s(i)).*sin((s(i) + 0.5)*theta) - ...
                    theta * sin(pi*s(i)).*cos((s(i) + 0.5)*theta);
    else
        xp(i) = pi*cos(pi*s(i)).*cos((s(i) - 0.5)*theta) - ...
                    theta * sin(pi*s(i)).*sin((s(i) - 0.5)*theta);
        yp(i) = pi*cos(pi*s(i)).*sin((s(i) - 0.5)*theta) + ...
                    theta * sin(pi*s(i)).*cos((s(i) - 0.5)*theta);
    end
end

zp = xp + 1i * yp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zpp = corner_zpp(s, theta)

xpp = zeros(size(s));
ypp = zeros(size(s));

for i = 1:length(s)
    if s(i) < 0
        xpp(i) = pi^2*sin(pi*s(i)).*cos((s(i)+0.5)*theta) + ...
                    pi*theta*cos(pi*s(i)).*sin((s(i)+0.5)*theta) + ...
                    theta*pi*cos(pi*s(i)).*sin((s(i)+0.5)*theta) + ...
                    theta^2*sin(pi*s(i)).*cos((s(i)+0.5)*theta);
        
        ypp(i) = pi^2*sin(pi*s(i)).*sin((s(i)+0.5)*theta) - ...
                    pi*theta*cos(pi*s(i)).*cos((s(i)+0.5)*theta) - ...
                    theta*pi*cos(pi*s(i)).*cos((s(i)+0.5)*theta) + ...
                    theta^2*sin(pi*s(i)).*sin((s(i)+0.5)*theta);
    else
        xpp(i) = -pi^2*sin(pi*s(i)).*cos((s(i)-0.5)*theta) - ...
                    pi*theta*cos(pi*s(i)).*sin((s(i)-0.5)*theta) - ...
                    theta*pi*cos(pi*s(i)).*sin((s(i)-0.5)*theta) - ...
                    theta^2*sin(pi*s(i)).*cos((s(i)-0.5)*theta);
        
        ypp(i) = -pi^2*sin(pi*s(i)).*sin((s(i)-0.5)*theta) + ...
                    pi*theta*cos(pi*s(i)).*cos((s(i)-0.5)*theta) + ...
                    theta*pi*cos(pi*s(i)).*cos((s(i)-0.5)*theta) - ...
                    theta^2*sin(pi*s(i)).*sin((s(i)-0.5)*theta);
    end
end

zpp = xpp + 1i * ypp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = square_z(t)

z = zeros(size(t));

for k = 1 : numel(t)
    if t(k) < pi / 2
        z(k) = 1i * t(k) / (pi / 2);
        
    else
        if t(k) < pi
            z(k) =  (t(k) - pi / 2) / (pi / 2) + 1i;
            
        else
            if t(k) < 3 * pi / 2
                z(k) = 1 - 1i * (pi - t(k)) / (pi / 2);
                
                
            else
                z(k) = 1  - (t(k) - 3 * pi / 2 )/ (pi / 2);
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zp = square_zp(t)

zp = zeros(size(t));

for k = 1 : numel(t)
    if t(k) < pi / 2
        zp(k) = -1i;
        
    else
        if t(k) < pi
            zp(k) = -1;
            
        else
            if t(k) < 3 * pi / 2
               zp(k) = 1i;                
                
            else
                zp(k) = 1;
            end
        end
    end
end

zp = zp / (pi / 2);



end
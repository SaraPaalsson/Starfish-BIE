function Disc = test_logkernelSQ()
% Function to test SQ for log-kernel and construct Lmod matrix.
%
% Integral to compute is I = Im( Int( f(tau)*log(tau-z) dtau ).

% Set parameters
Disc.Npanels = 4;
Disc.N = Disc.Npanels*16;
Disc.Param.interfparam = 'circle';

% Initialise discretisations
[Disc, ~] = init(Disc);

% figure(1); clf
% plot(Disc.z); hold on; axis equal;
% % plot(Dom.z+1i*eps,'*');
% quiver(real(Disc.z),imag(Disc.z),real(Disc.zp),imag(Disc.zp));

% Define function f
Disc.f = ones(Disc.N,1)*(1+1i)/10;

% Compute the integral I for z on boundary 
Disc.u  = compute_Ibndry(Disc);

figure(2); clf
plot(Disc.tpar,real(Disc.u)); 
hold on; 
plot(Disc.tpar,imag(Disc.u),'-.')

% Apply special quadrature for u on boundary
usq = compute_SQ(Disc);

plot(Disc.tpar,real(usq)); 
hold on; 
plot(Disc.tpar,imag(usq),'-.')


end


function L = compute_Lmod(zsrc,ztar)
% Compute SQ corrections for all points ztar from panel with points zsrc

Ntar = length(ztar);

pidx = (0:15);
V = ztar.^pidx;

for j=1:Ntar
    
    
end


end



function usq = compute_SQ(Disc)
% Computes onsurface SQ. This should be changed into 

load LmodFull.mat

wazp = Disc.w.*abs(Disc.zp);

usq = Disc.u;
usq = usq - Disc.f.*wazp.*log(pi/Disc.Npanels*abs(Disc.zp));

% Go through all panels. For each panel, modify the panel points and four
% points to each side of it.
for k=1:Disc.Npanels 
    idx = Disc.idx(1,1)+mod((k-1)*16+(-3:20)+Disc.Npanels*16-1,Disc.Npanels*16);
    idx2 = Disc.idx(1,1)-1+((k-1)*16+(1:16));
    
    Lmodtest = compute_Lmod(Disc.z(idx2),Disc.z(idx));
    
    usq(idx) = usq(idx) - Lmod*(Disc.f(idx2).*wazp(idx2));
end


end

function u = compute_Ibndry(Disc)
% Compute integral I = Im( Int( f(tau)*log(tau-z) dtau ) for z on bndry

u = zeros(Disc.N,1);
for j = 1:Disc.N
    tmp = Disc.f.*log(Disc.z-Disc.z(j)).*Disc.zp.*Disc.w; 
    tmp(j) = 0;
    u(j) = sum(tmp);
end

end


function [Disc, Dom] = init(Disc)
switch Disc.Param.interfparam
    case 'circle'
        r = 2;
        Disc.Param.tau = @(t) r*(cos(t) + 1i*sin(t));
        Disc.Param.taup = @(t) r*(-sin(t) + 1i*cos(t));
        Disc.Param.taupp = @(t) r*(-cos(t) - 1i*sin(t));
end
% Obtan GL-16 nodes and weights
[z, w,tpar,panels,pointPanel] = gl16(Disc.Npanels,Disc.Param.tau);
Disc.z = z; Disc.w = w; Disc.tpar = tpar; Disc.panels = panels;
Disc.pointPanel = pointPanel; 
Disc.idx = [1 Disc.N];
Disc.zp = Disc.Param.taup(Disc.tpar);
Disc.zpp = Disc.Param.taupp(Disc.tpar);

Dom.reld = linspace(0.25,0.99,20)';
Dom.z = Disc.z(end-5)*Dom.reld;

end


function [zDrops, wDrops, tpar,panels, pointPanel] = gl16(Npanels,tau)
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
    pointPanel((i-1)*16+1:16*i) = i;
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
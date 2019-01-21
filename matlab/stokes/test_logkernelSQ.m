function Disc = test_logkernelSQ()
% Function to test SQ for log-kernel and construct Lmod matrix.
%
% Integral to compute is I = Im( Int( f(tau)*log(tau-z) dtau ).

% Try to compute Lmod and compare with Rikard's:
L = compute_Lmod();

L

% % % Set parameters
% % Disc.Npanels = 4;
% % Disc.N = Disc.Npanels*16;
% % Disc.Param.interfparam = 'circle';
% % 
% % % Initialise discretisations
% % [Disc, ~] = init(Disc);
% % 
% % % figure(1); clf
% % % plot(Disc.z); hold on; axis equal;
% % % % plot(Dom.z+1i*eps,'*');
% % % quiver(real(Disc.z),imag(Disc.z),real(Disc.zp),imag(Disc.zp));
% % 
% % % Define function f
% % Disc.f = ones(Disc.N,1)*(1+1i)/10;
% % 
% % % Compute the integral I for z on boundary 
% % Disc.u  = compute_Ibndry(Disc);
% % % 
% % % figure(2); clf
% % % plot(Disc.tpar,real(Disc.u)); 
% % % hold on; 
% % % plot(Disc.tpar,imag(Disc.u),'-.')
% % 
% % % % Apply special quadrature for u on boundary
% % % usq = compute_SQ(Disc);
% % % 
% % % plot(Disc.tpar,real(usq)); 
% % % hold on; 
% % % plot(Disc.tpar,imag(usq),'-.')

end



function W = Wlogmat(tt)
npt = length(tt);
V = fliplr(vander(tt));
p = zeros(npt+1,1);
q = zeros(npt,1);
c = (1-(-1).^(1:npt))./(1:npt);
W = zeros(npt);
for j = 1:npt
    p(1) = log(abs((1-tt(j))/(1+tt(j))));
    for k=1:npt
        p(k+1) = tt(j)*p(k)+c(k);
    end
    q(1:2:npt-1) = log(abs(1-tt(j)^2))-p(2:2:npt);
    q(2:2:npt) = p(1)-p(3:2:npt +1);
    q = q./(1:npt)';
    wqj = (V')\q;
    W(j,:) = wqj';    
end
end


function L = compute_Lmod()
% Compute Lmod matrix for a straight panel of 16 G-L points with the panel
% points as targets

% Discretise panel
[zsrc,worig] = gaussleg(16,[-1 1]); % Create nodes and weights on the straight panel
% This panel goes between the points
tau1 = -1; 
tau2 = 1;
mid = (tau2+tau1)/2;
len = tau2-tau1;

W = Wlogmat(zsrc)

% Assign target points
[z1,~] = gaussleg(16,[-3 -1]); %left panel
[z2,~] = gaussleg(16,[1 3]); %right panel
% ztar = [z1(end-3:end); zsrc; z2(1:4)]

ztar = zsrc;

L = zeros(length(ztar),16);

for j=1:length(ztar) % Go through all points and compute wv
    % Assign target point
    z = ztar(j);
    nz = 2*(z-mid)/len;
         
    lg1 = log(abs(1-nz));
    lg2 = log(abs(1+nz));
    

    % Compute p0 analytically. NB p0 corresponds to p(1)
    p = zeros(17,1); 
    p(1) = lg1-lg2;
    
    % Compute p and r recursively
    r = zeros(16,1);
    gamma = len/2;
    for k = 2:17
        % Update p_k
        p(k) = nz*p(k-1) + (1-(-1)^(k-1))/(k-1);
        
        % Update r_{k-1}
        r(k-1) = 1/(k-1)*(lg1-(-1)^(k-1)*lg2-p(k));
        
    end
    
    wv = vandernewton(zsrc,r,16);
    wv = wv./worig;
    
    tmp = ~(zsrc==ztar(j));
    wv(tmp) = wv(tmp) - log(abs(zsrc(tmp)-ztar(j)));
    
    L(j,:) = wv.';
end


fprintf('Difference between L (mine) and W (AK): %e \n',norm(L-W,inf))

load('LmodFull.mat');
Lmod2 = Lmod(5:end-4,:);
fprintf('Difference between L (mine) and Lmod (Rikard): %e \n',norm(L-Lmod2,inf))

% disp('L:')
% L
% disp('');
% disp('Lmod2:')
% Lmod2
% disp('');

end

function b = vandernewton(T,b,n)

for k=1:n-1
    for i=n:-1:k+1
        b(i) = b(i) - T(k)*b(i-1);
    end
end

for k=n-1:-1:1
    for i=k+1:n
        b(i) = b(i)/(T(i)-T(i-k));
    end
    for i=k:n-1
        b(i) = b(i) - b(i+1);
    end
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
    
    Lmodtest = compute_Lmod(Disc.z(idx2),Disc.z(idx),Disc,k);
    
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
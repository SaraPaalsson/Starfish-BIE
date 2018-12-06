function [mu_fine, dom_fine] = ...
                recover_fine_density(dom_coarse, R, mu_tilde, star_indices, nsub)

mu_fine = zeros(length(mu_tilde) + 2 * nsub * 16, 2);
star_indices = star_indices(1 : end / 2);

b_star_indices = [17 : 80, 113 : 176];

mu_tilde = [real(mu_tilde), imag(mu_tilde)];

mu_tilde_c = [mu_tilde(star_indices,1); mu_tilde(star_indices,2)];

mu_fine_left1 = [];
mu_fine_right1 = [];
mu_fine_left2 = [];
mu_fine_right2 = [];

T = Tinit16;
W = Winit16;

Ib = eye(192);

[IP,~]=IPinit(T,W);
Pbc =blkdiag(eye(16),IP ,IP,eye(32), IP, IP, eye(16));

for i = nsub + 1 : -1 : 2
    
    [z, zp, zpp, jac] = geom_local_init(dom_coarse, T, W, nsub, i - 1);
    
    Kcirc = assemble_kernel(z, zp, zpp, jac) - Ib;
    Kcirc(b_star_indices, b_star_indices) = 0;
    
    FR = Ib + Kcirc;
    FR(b_star_indices, b_star_indices) = inv(R(:, :, i - 1));
    C = (Ib - Kcirc * inv(FR) ) *  Pbc;
    
    rho_vec_b = C * mu_tilde_c;
    
    mu_fine_left1 = [mu_fine_left1 ; rho_vec_b(1:16)];
    mu_fine_right1 = [rho_vec_b(81:96); mu_fine_right1 ];
    
    mu_fine_left2 = [mu_fine_left2 ; rho_vec_b(97 : 112)];
    mu_fine_right2 = [rho_vec_b(177 : 192); mu_fine_right2 ];
    
    mu_tilde_c = rho_vec_b(b_star_indices);
end

% middle two fine panels must be computed separetely
mu_fine_c = R(:,:,1) * mu_tilde_c;

mu_fine_left1 = [mu_fine_left1; mu_fine_c(1 : 32)];
mu_fine_right1 = [mu_fine_c(33 : 64); mu_fine_right1];

mu_fine_left2 = [mu_fine_left2; mu_fine_c(65 : 96)];
mu_fine_right2 = [mu_fine_c(97 : end); mu_fine_right2];

mu_fine_left = [mu_fine_left1, mu_fine_left2];
mu_fine_right = [mu_fine_right1, mu_fine_right2];

% place into correct indices
for i = 1 : 2
    
    starind_right = star_indices(end/2+1:end);   
    start_indice = 1;
    
    % read off density on Gamma^\circ
    mu_fine(start_indice : 16 * (dom_coarse.Npanels - 4) / 2, i) = ...
        mu_tilde(start_indice : 16 * (dom_coarse.Npanels - 4) / 2 ,i);
    
    start_indice =  16 * (dom_coarse.Npanels - 4) / 2 + 1;
    
    % add left corner
    mu_fine(start_indice : start_indice + size(mu_fine_right, 1) - 1, i) =...
        mu_fine_left(:,i);
    
    % add right corner
    start_indice = start_indice +  size(mu_fine_left, 1);
    mu_fine(start_indice : start_indice + size(mu_fine_right,1) - 1, i)...
        = mu_fine_right(:,i);
    
    % read off density on Gamma^\circ
    start_indice = start_indice +  size(mu_fine_right, 1);
    mu_fine(start_indice : end, i) = mu_tilde(starind_right(end) + 1: end, i);
    
end

mu_fine = mu_fine(:,1)+ 1i * mu_fine(:,2);

[z_fine, ~, jac_fine, s_fine, panel_breaks_s_fine] = ...
    geom_fine_init(dom_coarse, nsub, T, W);

% set up fine domain
ind_neg_coarse = 1 :  dom_coarse.N / 2 - 16;
ind_pos_coarse = dom_coarse.N / 2 + 17 : dom_coarse.N;
s_neg = dom_coarse.tpar(ind_neg_coarse);
s_pos = dom_coarse.tpar(ind_pos_coarse);
s_fine = [s_neg; s_fine; s_pos];
z_fine = [dom_coarse.tau(s_neg); z_fine; dom_coarse.tau(s_pos)];
jac_fine = [dom_coarse.wDrops(ind_neg_coarse); jac_fine; dom_coarse.wDrops(ind_pos_coarse)];

dom_fine = dom_coarse;
dom_fine.Npanels = dom_coarse.Npanels + 2 * nsub;
dom_fine.zDrops = z_fine;
dom_fine.tpar = s_fine;
dom_fine.panels = [dom_coarse.panels(1 : dom_coarse.Npanels / 2 - 1);panel_breaks_s_fine;...
                    dom_coarse.panels(dom_coarse.Npanels / 2 + 2 : end)];
dom_fine.wDrops = jac_fine;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = assemble_kernel(z,zp,zpp,W)

X = eye(2 * length(z));
K = zeros(size(X));

for i = 1 : 2 * length(z)
    K(:,i) = mubie_gmres(X(:,i),z,zp,zpp,W,[],[]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z, zp, jac, s, panel_breaks_s] = ...
    geom_fine_init(dom, nsub, T, W)

nq = length(W);
np = 2 * (nsub + 1) * nq;

s = zeros(np, 1);
z = zeros(np, 1);
zp = zeros(np, 1);
jac = zeros(np, 1);

coarse_panel_length = 1 / dom.Npanels;
panel_breaks_s = zeros(2 * nsub + 2,1);

for i = 1 : nsub
    indices_neg = (i - 1)*nq + 1 : i*nq; 
    indices_pos = np - i*nq + 1 : np - (i - 1)*nq;
    indices = [indices_neg, indices_pos];
    
    s_start = -coarse_panel_length / 2^(i - 1);
    s_end = -coarse_panel_length / 2^i;

    panel_breaks_s(i) = s_start;
    
    s(indices_neg) = T*(s_end - s_start)/2 + (s_start + s_end)/2;
    
    s_end = coarse_panel_length / 2^(i - 1);
    s_start = coarse_panel_length / 2^i;

    panel_breaks_s(2 * nsub - i + 3) = s_start;
   
    s(indices_pos) = T*(s_end - s_start)/2 + (s_start + s_end)/2;
    
    Wfine = [W; W] * (s_end - s_start)/2;
    
    z(indices) =  dom.tau(s(indices));
    zp(indices) = dom.taup(s(indices));
    
    jac(indices) = Wfine;
    
end

% the final fine panel
indices_neg = nsub*nq + 1 : (nsub + 1)*nq;
indices_pos = (nsub + 1)*nq + 1 : (nsub + 2)*nq;
indices = [indices_neg, indices_pos];

s_start = -coarse_panel_length / 2^nsub;

panel_breaks_s(nsub + 1) = s_start;

s_end = 0;

s(indices_neg) = T*(s_end - s_start)/2 + (s_start + s_end)/2;

s_end = coarse_panel_length / 2^nsub;
s_start = 0;

panel_breaks_s(nsub + 2) = s_start;

s(indices_pos) = T*(s_end - s_start)/2 + (s_start + s_end)/2;

Wfine = [W; W] * (s_end - s_start)/2;

z(indices) =  dom.tau(s(indices));
zp(indices) = dom.taup(s(indices));

jac(indices) = Wfine;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z, zp, zpp, jac] = ...
                geom_local_init(dom, T, W, nsub, level)
            
denom=2^(nsub-level)*dom.Npanels;

% construct s
s=[T/2 - 1.5; T/4 - 0.75; T/4 - 0.25; ...
    T/4 + 0.25; T/4 + 0.75; T/2 + 1.5]/denom;
W=[W/2; W/4; W/4; W/4; W/4; W/2]/denom;

z = dom.tau(s);
zp = dom.taup(s);
zpp = dom.taupp(s);

jac = W;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IP,IPW]=IPinit(T,W)
% *** create prolongation matrices ***
A=ones(16);
AA=ones(32,16);
T2=[T-1;T+1]/2;
W2=[W;W]/2;
for k=2:16
    A(:,k)=A(:,k-1).*T;
    AA(:,k)=AA(:,k-1).*T2;
end
IP=AA/A;
IPW=IP.*(W2*(1./W)');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=Tinit16
% *** 16-point Gauss-Legendre nodes ***
T=zeros(16,1);
T( 1)=-0.989400934991649932596154173450332627;
T( 2)=-0.944575023073232576077988415534608345;
T( 3)=-0.865631202387831743880467897712393132;
T( 4)=-0.755404408355003033895101194847442268;
T( 5)=-0.617876244402643748446671764048791019;
T( 6)=-0.458016777657227386342419442983577574;
T( 7)=-0.281603550779258913230460501460496106;
T( 8)=-0.095012509837637440185319335424958063;
T( 9)= 0.095012509837637440185319335424958063;
T(10)= 0.281603550779258913230460501460496106;
T(11)= 0.458016777657227386342419442983577574;
T(12)= 0.617876244402643748446671764048791019;
T(13)= 0.755404408355003033895101194847442268;
T(14)= 0.865631202387831743880467897712393132;
T(15)= 0.944575023073232576077988415534608345;
T(16)= 0.989400934991649932596154173450332627;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W=Winit16
% *** 16-point Gauss-Legendre weights ***
W=zeros(16,1);
W( 1)= 0.027152459411754094851780572456018104;
W( 2)= 0.062253523938647892862843836994377694;
W( 3)= 0.095158511682492784809925107602246226;
W( 4)= 0.124628971255533872052476282192016420;
W( 5)= 0.149595988816576732081501730547478549;
W( 6)= 0.169156519395002538189312079030359962;
W( 7)= 0.182603415044923588866763667969219939;
W( 8)= 0.189450610455068496285396723208283105;
W( 9)= 0.189450610455068496285396723208283105;
W(10)= 0.182603415044923588866763667969219939;
W(11)= 0.169156519395002538189312079030359962;
W(12)= 0.149595988816576732081501730547478549;
W(13)= 0.124628971255533872052476282192016420;
W(14)= 0.095158511682492784809925107602246226;
W(15)= 0.062253523938647892862843836994377694;
W(16)= 0.027152459411754094851780572456018104;

end

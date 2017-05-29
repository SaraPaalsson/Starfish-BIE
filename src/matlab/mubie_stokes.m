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

Afunc = @(x) mubie_gmres(x,zDrops,zpDrops,zppDrops,wDrops,[]);

w = gmres(Afunc,rhs,[],1e-13,100);
wr = w(1:N);
wi = w(N+1:end);
mu_stokes = wr + 1i*wi;


% vk = wDrops; tauk = zDrops; tau_primk = zpDrops; Ntau = N;
% 
% A = eye(2*N,2*N);
% rhs = zeros(2*N,1);
% G = @(x,xp,y) imag(xp./(x-y));
% F = @(x,xp,y) imag(xp.*(conj(x)-conj(y)))./((conj(x)-conj(y)).^2);
% for i=1:N
%     A(i,1:N) = 1/pi*(vk.*G(tauk,tau_primk,tauk(i)))- ...
%         1/pi*(vk.*real(F(tauk,tau_primk,tauk(i))));
%     A(i,Ntau+1:end) = -1/pi*(vk.*imag(F(tauk,tau_primk,tauk(i))));
%     A(i,i) = 1;
%     A(i,i+Ntau) = 0;
%     rhs(i) = imag(h(i));
%     
%     A(i+Ntau,1:Ntau) = -1/pi*(vk.*imag(F(tauk,tau_primk,tauk(i))));
%     A(i+Ntau,Ntau+1:end) = 1/pi*(vk.*G(tauk,tau_primk,tauk(i))) + ...
%         1/pi*(vk.*real(F(tauk,tau_primk,tauk(i))));
%     A(i+Ntau,i+Ntau) = 1;
%    A(i+Ntau,i) = 0;
%     rhs(i+Ntau) = -real(h(i));
% end
% w = gmres(A,-rhs,[],1e-10,Ntau);
% wr = w(1:Ntau);
% wi = w(Ntau+1:end);
% mu_stokes = wr + 1i*wi;

end
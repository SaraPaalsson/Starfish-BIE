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

end
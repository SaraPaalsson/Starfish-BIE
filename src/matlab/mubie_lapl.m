function mu_lapl = mubie_lapl(N,zDrops,zpDrops,zppDrops,wDrops,RHS)
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
mu_lapl = gmres(A,b,[],1e-13,N);
% mu_lapl = A\b;

end
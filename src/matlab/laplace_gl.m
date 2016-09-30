function [u,mu] = laplace_gl(tauk,tau_primk,tau_primprimk,zp,z,wk)

N = length(tauk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate source density mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = eye(N,N);
for i=1:N
    A(i,:) = A(i,:) + (1/pi*wk.*imag(tau_primk./(tauk-tauk(i))))';
end
d = 1 + 1/pi*wk.*imag(tau_primprimk./(2*tau_primk));
b = 2*imag(tauk.^2./(tauk-zp));
A(logical(eye(size(A)))) = 0;
A = A + diag(d,0);
%mu = gmres(A,b,[],1e-10,Ntau);
mu = A\b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute velocity u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u = zeros(size(z));
for j=1:length(z)
    for k=1:N
        if z(j) ~= tauk(k)
             u(j) = u(j) + 1/(2*pi)*wk(k)*mu(k)*imag(tau_primk(k)/(tauk(k)-z(j)));
        else
            u(j) = u(j) + 1/(2*pi)*wk(k)*mu(k)*imag(tau_primprimk(k)/(tau_primk(k)));
        end
    end
end

end
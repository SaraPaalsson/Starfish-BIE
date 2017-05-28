function u = laplace_trap(tauk,tau_primk,tau_primprimk,zp,z)

N = length(tauk);
% ---------------------------
% Set up trapezoidal weights
% % ---------------------------
dt = 2*pi/(N-1); 
wk = ones(N,1)*dt;
wk(1) = wk(1)/2;
wk(end) = wk(end)/2;

% ----------------------------
% Calculate A, b, mu from BIE
% ----------------------------
A = eye(N,N);
d = zeros(N,1);
b = zeros(N,1);
for i=1:N
    A(i,:) = A(i,:) + (1/pi*wk.*imag(tau_primk./(tauk-tauk(i))))';
end
d = 1 + 1/pi*wk.*imag(tau_primprimk./(2*tau_primk));
b = 2*imag(tauk.^2./(tauk-zp));
A(logical(eye(size(A)))) = 0;
A = A + diag(d,0);
A(1,N) = 1/pi*wk(end)*imag(tau_primprimk(1)/(2*tau_primk(1)));
A(N,1) = 1/pi*wk(1)*imag(tau_primprimk(1)/(2*tau_primk(1)));

mu = A\b;

% ----------------------------
% Calculate u from mu
% ----------------------------
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

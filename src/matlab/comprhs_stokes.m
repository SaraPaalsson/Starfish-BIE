function RHS = comprhs_stokes(x,x0,f)
% Compute bc for Stokes, using sum of stokeslets

% Number of point forces
Nf = size(f,1);

RHS = 0;
for j=1:Nf
    xhat = x-x0;
    r = abs(xhat);
    xh = real(xhat); yh = imag(xhat);
    f1 = real(f(j)); f2 = imag(f(j));
    S11 = -log(r) + (xh.^2)./(r.^2);
    S12 = (xh.*yh)./(r.^2);
    S21 = S12;
    S22 = -log(r) + (yh.^2)./(r.^2);
    RHS1 = S11*f1 + S12*f2;
    RHS2 = S21*f1 + S22*f2;
    RHS = RHS + (RHS1 + 1i*RHS2)/(4*pi);
end

end
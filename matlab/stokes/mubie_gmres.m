function Ax = mubie_gmres(x,z,zp,zpp,W,star_indices,Rhat)

if isempty(star_indices)
    Rhat = eye(length(x));
end

Rx = Rhat * x;
Nf = length(x)/2;

mu =  x(1:Nf) + 1i*x(Nf+1:end);
Rmu = Rx(1:Nf) + 1i*Rx(Nf+1:end);

M1_jj = imag(0.5*W.*zpp./zp);
M2_jj = imag(0.5*W.*zpp.*conj(zp))./(conj(zp).^2);

Aw = zeros(size(mu));
for j=1:length(mu)
    M1 = Rmu.*imag(W.*zp./(z-z(j)));
    M1(j) = Rmu(j)*M1_jj(j);
    
    M2 = conj(Rmu).*imag(W.*zp.*conj(z-z(j)))./(conj(z-z(j)).^2);
    M2(j) = conj(Rmu(j))*M2_jj(j);
    
    if ~isempty(star_indices) && max(j == star_indices)
        M1(star_indices(1:end/2)) = 0;
        M2(star_indices(1:end/2)) = 0;
    end
    
    Aw(j) = mu(j) + 1/pi*sum(M1) - 1/pi*sum(M2);
    
end

Ax = [real(Aw); imag(Aw)];

end
function Ax = mubie_gmres(x,z,zp,zpp,W)

Nf = length(x)/2;
mu = x(1:Nf) + 1i*x(Nf+1:end);


M1_jj = imag(0.5*W.*zpp./zp);
M2_jj = imag(0.5*W.*zpp.*conj(zp))./(conj(zp).^2);


Aw = zeros(size(mu));
for j=1:length(mu)
    M1 = mu.*imag(W.*zp./(z-z(j)));
    M1(j) = mu(j)*M1_jj(j);
    
    M2 = conj(mu).*imag(W.*zp.*conj(z-z(j)))./(conj(z-z(j)).^2);
    M2(j) = conj(mu(j))*M2_jj(j);
    
    Aw(j) = mu(j) + 1/pi*sum(M1) - 1/pi*sum(M2);
    
end


Ax = [real(Aw); imag(Aw)];


end
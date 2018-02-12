function u = compu_stokes(mu, z, zDrops, zpDrops, zppDrops, wDrops,evalinterf)
% Compute the velocity field by BIE for stokes equations, for all domain
% points z.
%
% OBS. We assume z not on boundary!
%
% Created 2017-05-28


if evalinterf
    z = zDrops;
    u = zeros(size(z));
    
    M1_jj = imag(0.5*wDrops.*zppDrops./zpDrops);
    M2_jj = 0.5*imag(wDrops.*zppDrops.*conj(zpDrops))./(conj(zpDrops).^2);
    
    for j=1:length(z)
        
        M1 = mu.*wDrops.*imag(zpDrops./(zDrops-z(j)));
        M1(j) = mu(j)*M1_jj(j);
        M2 = conj(mu).*wDrops.*imag(zpDrops.*conj(zDrops-z(j)))./(conj(zDrops-z(j)).^2);
        M2(j) = conj(mu(j))*M2_jj(j);
        
        u(j) = -1i/pi*sum(M1) + 1i/pi*sum(M2);
    end
    
    
else
    u = zeros(size(z));
    
    parfor j=1:length(z)
        %     M1 = sum(mu.*wDrops.*real(zpDrops./(zDrops-z(j))));
        %     M2 = sum(conj(mu).*wDrops.*imag(zpDrops.*conj(zDrops-z(j)))./(conj(zDrops-z(j)).^2));
        %     u(j) = -1/pi*M1 + 1i/pi*M2;
        
        M1 = sum(mu.*wDrops.*imag(zpDrops./(zDrops-z(j))));
        M2 = sum(conj(mu).*wDrops.*imag(zpDrops.*conj(zDrops-z(j)))./(conj(zDrops-z(j)).^2));
        
        u(j) = - 1i/pi*M1 + 1i/pi*M2;
    end
    
end

end
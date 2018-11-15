function errest = esterrors_stokes(zdom,zpanels,zDrops,zpDrops, ...
    Npanels,wDrops,mu)

gamma = @(t) t;
gammap = @(t) t.^0;

% qth derivative of k_n(z):
knq = @(n, q, z) (-(2*n+1)./sqrt(z.^2-1)).^q .* 2*pi./(z+sqrt(z.^2-1)).^(2*n+1);

errest = zeros(size(zdom));
for j=1:length(zdom)
    
    zj = zdom(j);
    
    % find closes panel edge
    [~,panel_idx] = min(abs(zj-zpanels));
    
    tmp1 = 0; tmp2 = 0; tmp3 = 0;
    for panidx = [-1 0] + panel_idx
        panidx2 = mod(panidx-1,Npanels) + 1;
        pidx = (panidx2-1)*16 + (1:16);
        
        z = zDrops(pidx); %Points on panel
        
        mid = (zpanels(panidx2+1)+zpanels(panidx2))/2;
        len = zpanels(panidx2+1) - zpanels(panidx2);
        
        z0 = 2*(zj-mid)/len;
        nz = 2*(z-mid)/len;
%         
%         za = zpanels(panidx2);
%         zb = zpanels(panidx2+1);
%        
%         d = min(abs(zj-z)); %Distance zj to closest point on panel
%         
%         
%         n1 = real(zb-za); n2 = imag(zb-za);
%         v1 = real(zj-za); v2 = imag(zj-za);
%         z0r = -1+2*(n1*v1+n2*v2)/(n1*n1+n2*n2);
%         L = sum(abs(wDrops(pidx).*zpDrops(pidx))); %Panel length
%         z0 = z0r + 1i*d*2/L; %Remapped z0.
%
        z0 = abs(real(z0)) + 1i*abs(imag(z0));

        mum1 = max(abs(1/(pi*1i)*mu(pidx).*gammap(nz)));
        tmp1 = tmp1 + abs(knq(16,0,z0)*mum1);

        ntz = -1i*zpDrops(pidx)./abs(zpDrops(pidx));
        mum2 = max(abs(1/(2*pi)*mu(pidx).*conj(ntz).^2.*gammap(nz)));
        tmp2 = tmp2 + abs(knq(16,0,z0)*mum2);
        
        mum3 = max(abs(1/(2*pi)*mu(pidx).*gammap(nz)));
        tmp3 = tmp3 + abs(mum3*knq(16,1,z0)*(conj(gamma(conj(z0)))-conj(gamma(z0)))/gammap(z0));
        
    end
    errest(j) = tmp1 + tmp2 + tmp3;
    
end

end
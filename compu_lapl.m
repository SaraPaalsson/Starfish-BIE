function u = compu_lapl(N, mu_lapl, z, zDrops, zpDrops, wDrops)

u = zeros(size(z));
% for j=1:length(z)
%     for k=1:N
%         if z(j) ~= zDrops(k)
%              u(j) = u(j) + 1/(2*pi)*wDrops(k)*mu_lapl(k)*imag(zpDrops(k)/(zDrops(k)-z(j)));
%         else
%             u(j) = u(j) + 1/(2*pi)*wDrops(k)*mu_lapl(k)*imag(zppDrops(k)/(zpDrops(k)));
%         end
%     end
% end

for j=1:length(z)
   u(j) = sum(mu_lapl.*wDrops.*imag(zpDrops./(zDrops-z(j)))); 
end
u = 1/(2*pi)*u;

end
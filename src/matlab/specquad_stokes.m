function [u_spec] = specquad_stokes(u, mudens, Npanels, tau_panels, zDrops, ...
    zpDrops, wDrops, z, IP1, IP2, W16, W32,u_known)

u_spec = u;
% Calculate mid points and lengths of all panels
mid2 = (tau_panels(2:end)+tau_panels(1:end-1))/2;
len2 = tau_panels(2:end)-tau_panels(1:end-1);


% %For debug
% z16 = zeros(size(z)); z32 = z16; zinter = z16;

deBUGmode = 0;
if deBUGmode
figure(1); clf;
plot(zDrops,'k-'); hold on; axis equal;
end

nbr_z = length(z);
for il = 1:nbr_z %Go through all points z
  
if deBUGmode
il 

figure(1);
plot(z(il),'*')
end

    for k=1:Npanels

if deBUGmode   
k
end       
        if abs(z(il)-mid2(k)) < abs(len2(k)) %Check if z too close to any panel
            tz = zeros(16,1); tzp = zeros(16,1); tW = zeros(16,1); tmu = zeros(16,1);
            nzpan = zeros(16,1);
            orig32 = zeros(32,1); tW32 = zeros(32,1);
            
            % % Calculate p0
            nz = 2*(z(il)-mid2(k))/len2(k); %map zk to nz (z0 in paper)
            lg1 = log(1-nz);
            lg2 = log(-1-nz);
            div1 = 1/(1+nz);
            div2 = 1/(1-nz);
            q32 = zeros(32,1);
            q32(1) = -div1 - div2;
            
            j = (1:16)';
            indj = (k-1)*16+j;
            
            tz = zDrops(indj);
            nzpan = 2*(tz-mid2(k))/len2(k);  %map all nodes on the panel to [-1,1]
            
if deBUGmode     
figure(3); clf;
plot(nzpan,'k.-','MarkerSize',20);
hold on
plot(nzpan,0*nzpan,'r-')


figure(3)
plot(nz,'*')
end
            %Check if the point nz is between the panel and the real axis
            if real(nz) > -1 && real(nz) < 1
                if imag(nz) > 0 %above real axis, check if enclosed by axis and panel
                    furthercheck = 0;
                    
                    if sum(imag(nzpan) > imag(nz)) ~= 0
                        furthercheck = 1;
                    end
                    
                    if furthercheck
                        %interpol. nzpan to poly and check value for
                        %real(nz)
                        tmpT = real(nzpan);
                        tmpb = imag(nzpan);
                        
                        p = vandernewtonT(tmpT,tmpb,16);
                        
                        kk = (0:15)';
                        test = sum(p.*real(nz).^kk);
                        
                        if test > imag(nz) %Correct value of integral
                            lg1 = lg1 - pi*1i; 
                            lg2 = lg2 + pi*1i;
%                             lg1 = lg1 + pi*1i;
%                             lg2 = lg2 - pi*1i;

                        end
                    end
                else if imag(nz) < 0 %below the real axis, check enclosed
                        furthercheck = 0;
                        
                        if sum(imag(nzpan) < imag(nz)) ~= 0
                            furthercheck = 1;
                        end
                        
                        if furthercheck
                            tmpT = real(nzpan);
                            tmpb = imag(nzpan);
                            
                            p = vandernewtonT(tmpT,tmpb,16);
                            
                            kk = (0:15)';
                            test = sum(p.*real(nz).^kk);
                            
                            if test < imag(nz) %Correct value of integral
                                lg1 = lg1 + pi*1i; 
                                lg2 = lg2 - pi*1i;
%                                 lg1 = lg1 - pi*1i; 
%                                 lg2 = lg2 + pi*1i;
                            end
                        end
                    end
                end
            end
            p32 = zeros(32,1);
            p32(1) = lg1-lg2;
            
            
            tzp = zpDrops(indj);
            tmu = mudens(indj);
            tW = wDrops(indj);
            % % Calculate old contribution to u from panel
            oldsum = sum(-1i/pi*tW.*tmu.*imag(tzp./(tz-z(il))) + ...
                1i/pi*tW.*conj(tmu).*imag(tzp.*conj(tz-z(il)))./(conj(tz-z(il)).^2));
            testsum = sum(tW.*tzp./(tz-z(il)));
if deBUGmode 
abs(p32(1)-testsum)
end
            if abs(p32(1)-testsum) > 1e-13 %Standard 16-GL not good enough!
                % % Interpolate to 32-point GL quadrature
                tmu32 = IPmultR(tmu,IP1,IP2);
                tz32 = IPmultR(tz,IP1,IP2);
                tzp32 = IPmultR(tzp,IP1,IP2);
                plen = tW(1)/W16(1);
                
                tW32 = W32.*plen;
                orig32 = tW32./(tz32-z(il));
                o32sum = sum(tzp32.*orig32);
if deBUGmode
abs(o32sum-p32(1))
end
                if abs(o32sum-p32(1)) < 1e-13 %32 GL suffices!
                    
                    
                    I1 = -1i/pi*tW32.*tmu32.*imag(tzp32./(tz32-z(il)));
                    I2 = 1i/pi*tW32.*conj(tmu32).*imag(tzp32.*conj(tz32-z(il)))./(conj(tz32-z(il)).^2);
                    newsum = sum(I1) + sum(I2);
                    
                    u_spec(il) = u_spec(il) + (newsum-oldsum);
                else %32 GL not enough, use interpolatory quadrature instead
                    % Use interpolatory quadrature
                    
            

                    nzpan32 = IPmultR(nzpan,IP1,IP2);
%                     gamma = 0.5*len2(k);
                    gamma = 2/len2(k);
                    
                    
                    signc = -1;
                    for j=1:31 %Calculate pk:s...
                        p32(j+1) = nz*p32(j) + (1-signc)/(j);
                        signc = -signc;
                        q32(j+1) = nz*q32(j) + p32(j);
                    end
                    
                    ntz = -1i*tzp32./abs(tzp32); %OK
                    f1 = tmu32.*conj(ntz).^2; %OK
                    f2 = tmu32.*conj(tz32-z(il)); %OK

                    
                    p32coeff = vandernewton(nzpan32,p32,32);
                    q32coeff = vandernewton(nzpan32,q32,32);
                                        
                    new1 = -1i/(pi) * (transpose(tmu32)*imag(p32coeff(1:end)));
                    new2 = conj(sum(f1.*p32coeff(1:end)));
                    new3 = conj(gamma*sum(f2.*q32coeff(1:end)));

%                     newsum = new1 + (0.5/pi)*new2 - (0.5/pi)*new3;
                    newsum = new1 + (-1)/(2*pi)*new2 + (-1)/(2*pi)*new3;

    
                    u_spec(il) = u_spec(il) + (newsum-oldsum);
                    
                end
                
                
            end
            
        end
if deBUGmode
u(il)
u_spec(il)
u_known(il)

abs(u_known(il)-u(il))

abs(u_known(il)-u_spec(il))
end
    end
    
end

end




function b = vandernewton(T,b,n)

for k=1:n-1
    for i=n:-1:k+1
        b(i) = b(i) - T(k)*b(i-1);
    end
end

for k=n-1:-1:1
    for i=k+1:n
        b(i) = b(i)/(T(i)-T(i-k));
    end
    for i=k:n-1
        b(i) = b(i) - b(i+1);
    end
end

%
% for k=1:n-1
%     for i=n:-1:k+1
%         b(i) = b(i) - T(k)*b(i-1);
%     end
% end
% for k=n-1:-1:1
%     for i=k+1:n
%         b(i) = b(i)/(x(i)-x(i-k));
%     end
%     for i=k:n-1
%         b(i) = b(i) - b(i+1);
%
%     end
% end
%

end
% function b = vandernewton(x,b,n)
% for k=1:n-1
%     for i=n:-1:k+1
%         b(i) = b(i) - x(k)*b(i-1);
%     end
% end
% for k=n-1:-1:1
%     for i=k+1:n
%         b(i) = b(i)/(x(i)-x(i-k));
%     end
%     for i=k:n-1
%         b(i) = b(i) - b(i+1);
%
%     end
% end
% end




function [a] = vandernewtonT(T,b,n)
x = T;
c = b;
for k=1:n-1
    for i=n:-1:k+1
        c(i) = (c(i)-c(i-1))/(x(i)-x(i-k));
    end
end
a = c;
for k=n-1:-1:1
    for i=k:n-1
        a(i) = a(i)-x(k)*a(i+1);
    end
end
end

% function [a] = vandernewtonT(T,b,n)
% x = T;
% c = b;
% for k=1:n-1
%     for i=n:-1:k+1
%         c(i) = (c(i)-c(i-1))/(x(i)-x(i-k));
%     end
% end
% a = c;
% for k=n-1:-1:1
%     for i=k:n-1
%         a(i) = a(i)-x(k)*a(i+1);
%     end
% end
% end


% function b = vandernewton(T,b,n)
%
% for k = 2:n
%    for j=n:-1:k
%       b(j) = b(j) - T(k-1)*b(j-1);
%    end
% end
% for k=n:-1:1
%    for j=k+1:n
%       b(j) = b(j)/(T(j)-T(j-k));
%       b(j-1) = b(j-1) - b(j);
%    end
% end
% % for k = 2:n
% %    for j=n:-1:k+1
% %       b(j) = b(j) - T(k-1)*b(j-1);
% %    end
% % end
% % for k=n:-1:2
% %    for j=k+1:n
% %       b(j) = b(j)/(T(j)-T(j-k));
% %       b(j-1) = b(j-1) - b(j);
% %    end
% % end
% end

%
% function [T,b] = vandernewtonT(T,b,n)
% x = T;
% c = b;
% for k=1:n-1
%     for i=n:-1:k+1
%         c(i) = (c(i)-c(i-1))/(x(i)-x(i-k));
%     end
% end
% a = c;
% for k=n-1:-1:1
%     for i=k:n-1
%         a(i) = a(i)-x(k)*a(i+1);
%     end
% end
% %FEL FEL FEL
% % for k=1:n-1
% %    for j=n:-1:(k+2)
% %        b(j) = (b(j)-b(j-1))/(T(j)-T(j-k-1));
%    end
% end
% for k=n:-1:1
%     for j=k:n-1
%         b(j) = b(j) - T(k)*b(j+1);
%     end
% end
% A = zeros(n,n+1);
% A(:,end) = 1;
% k = n:-1:1; k = k(:);
% for j=1:n
%     A(:,j) = T.^(k(j));
% end
% % A(:,1:n-1) = T.^k;
% b = A\b;
% end

function [f32] = IPmultR(f16,IP1,IP2)
f32 = zeros(32,1);
for i=1:16
   t1 = 0;
   t2 = 0;
   ptr = i;
   for j=1:8
      t1 = t1 + IP1(ptr)*(f16(j)+f16(17-j)); %17
      t2 = t2 + IP2(ptr)*(f16(j)-f16(17-j)); %17
      ptr = ptr + 16;
   end
   f32(i) = t1+t2;
   f32(33-i) = t1-t2; %33
end
end
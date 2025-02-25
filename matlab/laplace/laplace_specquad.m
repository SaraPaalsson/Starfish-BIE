function [u_spec] = laplace_specquad(u, mudens, Npanels, tau_panels, zDrops, ...
    zpDrops, wDrops, z, IP1, IP2, W16, W32)

u_spec = u;
% Calculate mid points and lengths of all panels
mid2 = (tau_panels(2:end)+tau_panels(1:end-1))/2;
len2 = tau_panels(2:end)-tau_panels(1:end-1);


% %For debug
% z16 = zeros(size(z)); z32 = z16; zinter = z16;

nbr_z = length(z);
parfor i = 1:nbr_z %Go through all points z   
        for k=1:Npanels
             
%             if 1
            if abs(z(i)-mid2(k)) < abs(len2(k)) %Check if z too close to any panel
                tz = zeros(16,1); tzp = zeros(16,1); tW = zeros(16,1); tmu = zeros(16,1);
                nzpan = zeros(16,1); 
                orig32 = zeros(32,1); tW32 = zeros(32,1);
                
                % % Calculate p0
                nz = 2*(z(i)-mid2(k))/len2(k); %map zk to nz (z0 in paper)
                oldsum = 0; testsum = 0;
                lg1 = log(1-nz);
                lg2 = log(-1-nz);
                
                j = (1:16)';
                indj = (k-1)*16+j;
                tz = zDrops(indj);
                nzpan = 2*(tz-mid2(k))/len2(k);  %map all nodes on the panel to [-1,1]
                
                
                %Check if the point nz is between the panel and the real axis
                if real(nz) > -1 && real(nz) < 1
                    if imag(nz) > 0 %above real axis, check if enclosed by axis and panel
                        furthercheck = 0;
                        for j=1:16 %go through all nodes on panel
                            if imag(nzpan(j)) > imag(nz)
                                    furthercheck = 1;
                                    break;
                            end
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
                            end
                        end
                    else if imag(nz) < 0 %below the real axis, check enclosed
                            furthercheck = 0;
                            for j=1:16 %go through all nodes on panel
                                if imag(nzpan(j)) < imag(nz)
                                    furthercheck = 1;
                                    break;
                                end
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
                                end
                            end
                        end
                    end
                end
                p32 = zeros(32,1);
                p32(1) = lg1-lg2;
                % % Calculate old contribution to u from panel
                
                
                tzp = zpDrops(indj);
                tmu = mudens(indj);
                tW = wDrops(indj);
                oldsum = 1/(2*pi)*sum(tW.*tmu.*imag(tzp./(tz-z(i))));
                testsum = sum(tW.*tzp./(tz-z(i)));
      

%                 if 1
               if abs(p32(1)-testsum) > 1e-13 %Standard 16-GL not good enough!
                    % % Interpolate to 32-point GL quadrature
                    tmu32 = IPmultR(tmu,IP1,IP2);
                    tz32 = IPmultR(tz,IP1,IP2);
                    tzp32 = IPmultR(tzp,IP1,IP2);
                    plen = tW(1)/W16(1);

                    tW32 = W32.*plen;
                    orig32 = tW32./(tz32-z(i));
                    o32sum = sum(tzp32.*orig32);
                    
%                     if 0
                    if abs(o32sum-p32(1)) < 1e-13 %32 GL suffices!
                        
                        newsum = 1/(2*pi)*sum(tW32.*tmu32.*imag(tzp32./(tz32-z(i))));
                        
                        u_spec(i) = u_spec(i) + (newsum-oldsum);
                    else %32 GL not enough, use interpolatory quadrature instead
                        % Use interpolatory quadrature
                        
                        nzpan32 = IPmultR(nzpan,IP1,IP2);
                        signc = -1;
                        for j=1:31 %Calculate pk:s...
                            p32(j+1) = nz*p32(j) + (1-signc)/(j);%(1-signc)/j;
                            signc = -signc;
                        end
%                         c32 = vandernewtonT(nzpan32,tmu32,32);
%                         newsum = 0;
%                         for j=1:32
%                             newsum = newsum + 1/(2*pi) * imag(p32(j)*c32(j));
%                         end
%                         
                        p32coeff = vandernewton(nzpan32,p32,32);  
                        newsum2 = 1/(2*pi)*sum(imag(p32coeff.*tmu32));
                        
                        u_spec(i) = u_spec(i) + (newsum2-oldsum);
                                              
                    end
                   
                
                end
                
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
end


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
function [u_spec] = specquad_lapl(u, mudens, Npanels, panels, zDrops, ...
    zpDrops, wDrops, z, IP1, IP2, W16, W32)


u_spec = u;


% Calculate mid points and lengths of all panels
% mid2 = (tau(panels(2:end))+tau(panels(1:end-1)))/2;
% len2 = tau(panels(2:end))-tau(panels(1:end-1));
mid2 = (panels(2:end)+panels(1:end-1))/2;
len2 = panels(2:end)-panels(1:end-1);



nbr_z = length(z);
for i = 1:nbr_z %Go through all points z   
% for i=1:400
        for k=1:Npanels
            
            if abs(z(i)-mid2(k)) < abs(len2(k)) %Check if z too close to any panel

                % % Calculate p0
                nz = 2*(z(i)-mid2(k))/len2(k); %map zk to nz (z0 in paper)
                lg1 = log(1-nz);
                lg2 = log(-1-nz);
                
                tz = zDrops((k-1)*16+(1:16));
                nzpan = 2*(tz-mid2(k))./len2(k);

                %Check if the point nz is between the panel and the real axis
                if real(nz) > -1 && real(nz) < 1
                    if imag(nz) > 0 %above real axis, check if enclosed by axis and panel
                        furthercheck = sum(imag(nzpan) > imag(nz));

                        if furthercheck > 0
                            %interpol. nzpan to poly and check value for
                            %real(nz)
                            tmpT = real(nzpan);
                            tmpb = imag(nzpan);

                            p = vandernewtonT(tmpT,tmpb,16);
                            test = sum(p(1:16).*real(nz).^(0:15)');
                                
                                
                            if test > imag(nz) %Correct value of integral
%                                 lg1 = lg1 + pi*1i;
%                                 lg2 = lg2 - pi*1i;
                                lg1 = lg1 - pi*1i; %HMMM???
                                lg2 = lg2 + pi*1i;
                            end
                        end
                    else if imag(nz) < 0 %below the real axis, check enclosed

                            furthercheck = sum(imag(nzpan)<imag(nz));
                            
                            if furthercheck > 1
                                tmpT = real(nzpan);
                                tmpb = imag(nzpan);
  
                                p = vandernewtonT(tmpT,tmpb,16);
                                
                                test = sum(p(1:16).*real(nz).^(0:15)');


                                if test < imag(nz) %Correct value of integral
%                                     lg1 = lg1 - pi*1i;
%                                     lg2 = lg2 + pi*1i;
                                    lg1 = lg1 + pi*1i; %HMMM??
                                    lg2 = lg2 - pi*1i;
                                end
                            end
                        end
                    end
                end
                p32 = zeros(32,1);
                p32(1) = lg1-lg2;
               
                
                % % Calculate old contribution to u from panel
                tzp = zpDrops((k-1)*16+(1:16));
                tmu = mudens((k-1)*16+(1:16));
                tW = wDrops((k-1)*16+(1:16));
                oldsum = sum(tW.*tmu.*imag(tzp./(tz-z(i))))/2/pi;
                testsum = sum(tW.*tzp./(tz-z(i)));
                

                if abs(p32(1)-testsum) > 1e-13 %Standard 16-GL not good enough!
                    % % Interpolate to 32-point GL quadrature
                    tmu32 = IPmultR(tmu,IP1,IP2);
                    tz32 = IPmultR(tz,IP1,IP2);
                    tzp32 = IPmultR(tzp,IP1,IP2);
                    plen = tW(1)/W16(1);
                    tW32 = W32*plen;
                    orig32 = tW32./(tz32-z(i));
                    o32sum = sum(tzp32.*orig32);
                    
                    if abs(o32sum-p32(1)) < 1e-13 %32 GL suffices!
                        newsum = sum(tW32.*tmu32.*imag(tzp32./(tz32-z(i))))/2/pi;
%                         
                        u_spec(i) = u_spec(i) + (newsum-oldsum);
                        
                     
                        
                    else %32 GL not enough, use interpolatory quadrature instead
                        % Use interpolatory quadrature
                        
                        nzpan32 = IPmultR(nzpan,IP1,IP2);
                        signc = -1;
                        for j=1:31 %Calculate pk:s...
                            p32(j+1) = nz*p32(j) + (1-signc)/(j);%(1-signc)/j;
                            signc = -signc;
                        end

                        p32coeff = vandernewton(nzpan32,p32,32);
                        newsum2 = sum(imag(p32coeff.*tmu32))/2/pi;
                        
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
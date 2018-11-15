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

clc;
clf;
for i=-100:0.1:100
   %val = 1.0 + i+i*i*0.5 + i*i*i*(1.0/6.0) + i*i*i*i*(1.0/24.0);
   val = 1-i;
   if(abs(val))<=1.0
      plot(i,'*');
      hold on;
   end
end
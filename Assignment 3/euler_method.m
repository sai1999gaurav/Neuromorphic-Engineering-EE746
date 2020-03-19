function [V_t] = euler_method(Iapp,del_t)
global t1 t2 t3 ;
C = 300e-12;
gl = 30e-9;
El = -70e-3;
V_thresh = -20e-3;
R_p = 2e-3;
t_size = 2000;
V_t = zeros(3,t_size);
V_t(:,1) = El;
counter1 = 0;
counter2 = 0;
counter3 = 0;
for i = 2:t_size 
   V_t(1,i) = V_t(1,i-1) + (del_t/C)*(-gl*(V_t(1,i-1) - El) + Iapp(1,i-1));
   V_t(2,i) = V_t(2,i-1) + (del_t/C)*(-gl*(V_t(2,i-1) - El) + Iapp(2,i-1)); 
   V_t(3,i) = V_t(3,i-1) + (del_t/C)*(-gl*(V_t(3,i-1) - El) + Iapp(3,i-1));   
   counter1 = counter1 - 1;
   counter2 = counter2 - 1;
   counter3 = counter3 - 1;
   if(counter1 > 0)
       V_t(1,i) = El;
   end
   if(counter2 > 0)
       V_t(2,i) = El;
   end
   if(counter3 > 0)
       V_t(3,i) = El;
   end
   if (V_t(1,i) >=  V_thresh)
    V_t(1,i) = El;
    t1 = i;    
    counter1=20;
    sprintf('Neuron 1: %f ', i);
   end
   if (V_t(2,i) >=  V_thresh)
    V_t(2,i) = El;
    t2 = i;
    counter2=20;
    sprintf('Neuron 2: %f ', i);
   end
   if (V_t(3,i) >=  V_thresh)
    V_t(3,i) = El;
    t3 = i;
    counter3=20;
    sprintf('Neuron 3: %f ', i);
   end   
end
end
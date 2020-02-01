function [V_t] = euler_method(type,Iapp,del_t)
C = (100*(type - 2)*(type - 3) - 130*(type-1)*(type-3) + 100*(type-1)*(type-2))*1e-12;
gl = (5*(type - 2)*(type - 3) - 18*(type-1)*(type-3) + 5*(type-1)*(type-2))*1e-9;
El = (-35*(type - 2)*(type - 3) + 58*(type-1)*(type-3) - 29*(type-1)*(type-2))*1e-3;
V_thresh = -50e-3;
delta_volt = 2e-3;
a = (1*(type - 2)*(type - 3) - 4*(type-1)*(type-3) + (type-1)*(type-2))*1e-9;
tau_w = (15*(type - 2)*(type - 3) - 150*(type-1)*(type-3) + 60*(type-1)*(type-2))*1e-3;
b = (0*(type - 2)*(type - 3) - 120*(type-1)*(type-3) + 50*(type-1)*(type-2))*1e-12;
V_r = (-29*(type - 2)*(type - 3) + 50*(type-1)*(type-3) - 23*(type-1)*(type-2))*1e-3;
m = 5000;
U_t = zeros(m,1);
V_t = zeros(m,1);
V_t(1) = El;
for i = 2:m
   V_t(i) = V_t(i-1) + (del_t/C)*(-gl*(V_t(i-1) - El) + gl*delta_volt*exp((V_t(i-1) - V_thresh)/delta_volt) - U_t(i-1) + Iapp);
   U_t(i) = U_t(i-1) + (del_t/tau_w)*(a*(V_t(i-1) - El) -  U_t(i-1));
   if (V_t(i) >=  0)
    V_t(i) = V_r;
    U_t(i) = U_t(i) + b;
   end
   
end
end
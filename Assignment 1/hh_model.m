function [V_t, I_na, I_k, I_l] = hh_model(I_t, del_t)
s = 5000;
V_t = zeros(s,1);
I_na = zeros(s,1);
I_k = zeros(s,1);
I_l = zeros(s,1);
V_t(1) = -65; %Resting Potential
alpha_m = 0.1*(V_t(1) + 40)/(1 -  exp(-(V_t(1) + 40)/10));
alpha_n = 0.01*(V_t(1) + 55)/(1 -  exp(-(V_t(1) + 55)/10));
alpha_h = 0.07*exp(-0.05*(V_t(1) + 65));
beta_m = 4*exp(-0.0556*(V_t(1) + 65));
beta_n = 0.125*exp(-1*(V_t(1) + 65)/80);
beta_h = 1/(1+exp(-0.1*(V_t(1) + 35)));
V_t(1) = -65e-3;
%Initial Conditions
m_old = alpha_m/(alpha_m + beta_m); 
n_old = alpha_n/(alpha_n + beta_n);
h_old = alpha_h/(alpha_h + beta_h);
g_na = 120e-3;
g_k = 36e-3;
g_l = 0.3e-3;
E_na = 50e-3;
E_k = -77e-3;
E_l = -55e-3;
for i = 2:s
   m_new = m_func(V_t(i-1)*1e3, m_old, del_t);
   h_new = h_func(V_t(i-1)*1e3, h_old, del_t);
   n_new = n_func(V_t(i-1)*1e3, n_old, del_t);
   m_old = m_new;
   h_old = h_new;
   n_old = n_new;
   I_na(i) = g_na*m_old*m_old*m_old*h_old*(V_t(i-1) - E_na);
   I_k(i) = g_k*n_old*n_old*n_old*n_old*(V_t(i-1) - E_k);
   I_l(i) = g_l*(V_t(i-1) - E_l); 
   V_t(i) = volt_t(V_t(i-1), I_na(i-1), I_k(i-1), I_l(i-1),I_t(i-1),del_t); 
end

end

function [V_t] = volt_t(V_t_i, I_na_i, I_k_i, I_l_i, I_t_i, del_t)
C = 1e-6;
V_t = V_t_i + (del_t/C)*(I_t_i - I_na_i - I_k_i -I_l_i);
end

function [m_new] = m_func(V_t_i, m_old, del_t)
alpha = 0.1*(V_t_i + 40)/(1 -  exp(-(V_t_i + 40)/10));
beta = 4*exp(-0.0556*(V_t_i + 65));
m_new = m_old + del_t*(alpha*(1-m_old) - beta*m_old);
end

function [h_new] = h_func(V_t_i, h_old, del_t)
alpha = 0.07*exp(-0.05*(V_t_i + 65));
beta = 1/(1+exp(-0.1*(V_t_i + 35)));
h_new = h_old + del_t*(alpha*(1-h_old) - beta*h_old);
end

function [n_new] = n_func(V_t_i, n_old, del_t)
alpha = 0.01*(V_t_i + 55)/(1 -  exp(-(V_t_i + 55)/10));
beta = 0.125*exp(-1*(V_t_i + 65)/80);
n_new = n_old + del_t*(alpha*(1-n_old) - beta*n_old);
end


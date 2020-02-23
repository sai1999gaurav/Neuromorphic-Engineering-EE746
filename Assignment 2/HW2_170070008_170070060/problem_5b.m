function [syn_strength] = problem_5b(m_syn_stimulus_2, V_t_2)
Ns = 100;
w0 = 200;
sigma_w = 20;
T = 0.5;
dt = 1e-4;
l_t = 1;
lambda = l_t*dt;
time = (dt:dt:T);
% 5(b)
[syn_strength, num_iter] = problem_3(m_syn_stimulus_2,V_t_2, 50, 5);
fprintf('Number of Iterations: %d\n', num_iter);
fprintf('Weights: \n');
disp(syn_strength); % w1
end
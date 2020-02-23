%% 
%Problem 1
problem_1();
%% 
%Problem 2a
[m_syn_stimulus, V_t] = problem_2a();
%PART (b) -- SAME STIMULUS
V_t_2 = problem_2b(m_syn_stimulus);
% Problem 3 -- SAME STIMULUS
w0 = 50;    
sigma_w = 5;
[syn_strength, num_iter] = problem_3(m_syn_stimulus, V_t, w0, sigma_w);
fprintf('Number of Iterations: %d\n', num_iter);
fprintf('Weights: \n');
disp(syn_strength);
%PROBLEM 4
w0 = 250;
sigma_w = 25;
[syn_strength, num_iter] = problem_4(m_syn_stimulus,V_t_2, w0, sigma_w);
fprintf('Number of Iterations: %d\n', num_iter);
fprintf('Weights: \n');
disp(syn_strength);
%%
%Problem 5a
[m_syn_stimulus_1,V_t_1,m_syn_stimulus_2, V_t_2] = problem_5a();
%Problem 5b
syn_strength = problem_5b(m_syn_stimulus_2, V_t_2);
%Problem 5c
problem_5c(syn_strength, m_syn_stimulus_1, m_syn_stimulus_2, V_t_2, 1);
%Problem 5d
fprintf('5d starts');
problem_5d(m_syn_stimulus_1, V_t_1, m_syn_stimulus_2, V_t_2);
fprintf('5d end');
function problem_5c(syn_strength, m_syn_stimulus_1, m_syn_stimulus_2, V_t_2, numb)
Ns = 100;
w0 = 200;
sigma_w = 20;
T = 0.5;
dt = 1e-4;
l_t = 1;
lambda = l_t*dt;
time = (dt:dt:T);
t_size = length(time);
w1 = mean(syn_strength);
sigma_w1 = std2(syn_strength);
I0 = 1e-12;
tau = 15e-3;
tau_s = tau/4;
syn_strength_s1  = zeros(Ns,1);
Iapplied_1 = zeros(Ns,t_size);
for k = 1: Ns
    for i = 1:t_size
        for t = 1:i
            if (m_syn_stimulus_1(k) == t)
                Iapplied_1(k,i) = Iapplied_1(k,i) + I0*syn_strength(k,1)*(exp(-(i - t)*dt/tau) - exp(-(i - t)*dt/tau_s));  
            end
        end
    end
end 
I_total_1 = sum(Iapplied_1, 1); % Sum of all synapses 1 x 5000
%Membrane Potential (AEF RS) :HW 1 Q3
V_t_1 = euler_method(1,I_total_1,dt); %Same function as used in HW1
if (max(V_t_1) < -0.04) % No spike response for S1, then
 while (max(V_t_1) < -0.04)
 [syn_strength_s1, num_iter_s1] = problem_3(m_syn_stimulus_1, V_t_1, w1, sigma_w1);
 w2 = mean(syn_strength_s1);
 sigma_w2 = std2(syn_strength_s1);
 [syn_strength_s2, num_iter_s2] = problem_4(m_syn_stimulus_2, V_t_2, w2, sigma_w2);
 w1 = mean(syn_strength_s2);
 sigma_w1 = std2(syn_strength_s2);
 Iapplied_1 = zeros(Ns,t_size);
 for k = 1: Ns
    for i = 1:t_size
        for t = 1:i
            if (m_syn_stimulus_1(k) == t)
                Iapplied_1(k,i) = Iapplied_1(k,i) + I0*syn_strength(k,1)*(exp(-(i - t)*dt/tau) - exp(-(i - t)*dt/tau_s));  
            end
        end
    end
 end 
 I_total_1 = sum(Iapplied_1, 1); % Sum of all synapses 1 x 5000
 %Membrane Potential (AEF RS) :HW 1 Q3
 V_t_1 = euler_method(1,I_total_1,dt); %Same function as used in HW1

 end
end
if (numb == 1)
figure(9);
subplot(2,1,1);
plot(time.*1000, V_t_1, 'linewidth', 2);
title(sprintf('Neuron: RS Response Stimulus 1 Updated weights'));
xlabel('Time (in ms)');
ylabel('Spikes');
subplot(2,1,2);
plot(time.*1000, V_t_2, 'linewidth', 2);
title(sprintf('Neuron: RS Response Stimulus other one '));
xlabel('Time (in ms)');
ylabel('Spikes');
%fprintf('Number of Iterations: %d\n', num_iter);
fprintf('Weights: \n Mean: ');
disp(mean(syn_strength_s1)); % w1
fprintf('SD: %f', std2(syn_strength_s1) );

else % for 5d
figure(10);
subplot(2,1,1);
plot(time.*1000, V_t_1, 'linewidth', 2);
title(sprintf('Neuron: RS Response Stimulus 2 Updated weights'));
xlabel('Time (in ms)');
ylabel('Spikes');
subplot(2,1,2);
plot(time.*1000, V_t_2, 'linewidth', 2);
title(sprintf('Neuron: RS Response Stimulus other one '));
xlabel('Time (in ms)');
ylabel('Spikes');
%fprintf('Number of Iterations: %d\n', num_iter);
fprintf('Weights: \n Mean: ');
disp(mean(syn_strength_s1)); % w1
fprintf('SD: %f', std2(syn_strength_s1) );
   
end
    
end
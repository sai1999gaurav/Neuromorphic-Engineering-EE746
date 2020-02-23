function [m_syn_stimulus_1, V_t_1, m_syn_stimulus_2, V_t_2] = problem_5a()
%5(a)
Ns = 100;
w0 = 200;
sigma_w = 20;
T = 0.5;
dt = 1e-4;
l_t = 1;
lambda = l_t*dt;
syn_strength_old = w0 + sigma_w*randn(Ns,1); %Strength gaussian -- randn -- standard normal
time = (dt:dt:T);
poisson_output = makedist('Poisson', 'lambda', lambda);
t_size = length(time);
m_syn_stimulus_1 = zeros(Ns, 1);
m_syn_stimulus_2 = zeros(Ns, 1);
for i = 1:Ns % Each Synapse
    for j = 1:t_size % Each time
        spike = random(poisson_output,1);
        if (spike >= 1)
         m_syn_stimulus_1(i) = j; %spike-train
        end
        spike = random(poisson_output,1);
        if (spike >= 1)
         m_syn_stimulus_2(i) = j; %spike-train
        end
    end
end
I0 = 1e-12;
tau = 15e-3;
tau_s = tau/4;
Iapplied_1 = zeros(Ns,t_size);
Iapplied_2 = zeros(Ns, t_size);
for k = 1: Ns
    for i = 1:t_size
        for t = 1:i
            if (m_syn_stimulus_1(k) == t)
                Iapplied_1(k,i) = Iapplied_1(k,i) + I0*syn_strength_old(k,1)*(exp(-(i - t)*dt/tau) - exp(-(i - t)*dt/tau_s));  
            end
            if (m_syn_stimulus_2(k) == t)
                Iapplied_2(k,i) = Iapplied_2(k,i) + I0*syn_strength_old(k,1)*(exp(-(i - t)*dt/tau) - exp(-(i - t)*dt/tau_s));  
            end
        end
    end
end 
I_total_1 = sum(Iapplied_1, 1); % Sum of all synapses 1 x 5000
I_total_2 = sum(Iapplied_2, 1); % Sum of all synapses 1 x 5000
%Membrane Potential (AEF RS) :HW 1 Q3
V_t_1 = euler_method(1,I_total_1,dt); %Same function as used in HW1
V_t_2 = euler_method(1,I_total_2,dt); %Same function as used in HW1

figure(8);
subplot(2,1,1);
plot(time.*1000, V_t_1, 'linewidth', 2);
title(sprintf('Neuron: RS Response Stimulus 1'));
xlabel('Time (in ms)');
ylabel('Spikes');
subplot(2,1,2);
plot(time.*1000, V_t_2, 'linewidth', 2);
title(sprintf('Neuron: RS Response Stimulus 2'));
xlabel('Time (in ms)');
ylabel('Spikes');

end
function [V_t] = problem_2b(m_syn_stimulus)
w0 = 250;
sigma_w = 25;
Ns = 100;
T = 0.5;
dt = 1e-4;
l_t = 1;
lambda = l_t*dt;
syn_strength = w0 + sigma_w*randn(Ns,1); %Strength gaussian -- randn -- standard normal
time = (dt:dt:T);
t_size = length(time);
I0 = 1e-12;
tau = 15e-3;
tau_s = tau/4;
Iapplied = zeros(Ns,t_size);
for k = 1: Ns
    for i = 1:t_size
        for t = 1:i
            if (m_syn_stimulus(k) == t)
            Iapplied(k,i) = Iapplied(k,i) + I0*syn_strength(k,1)*(exp(-(i - t)*dt/tau) - exp(-(i - t)*dt/tau_s));  
            end
        end
    end
end 
I_total = sum(Iapplied, 1); % Sum of all synapses 1 x 5000
figure(3);
sgtitle(sprintf('Neuron Response with w0 = 250 and sigma_w = 25 (2(b))'));
subplot(2,1,1);
plot(time.*1000,I_total(1,:), 'Linewidth',2);
title('Synaptic Current');
ylabel('Iapplied');
xlabel('Time (in ms)');
%Membrane Potential (AEF RS) :HW 1 Q3
V_t = euler_method(1,I_total,dt); %Same function as used in HW1
subplot(2,1,2);
plot(time.*1000, V_t, 'linewidth', 2);
title(sprintf('Neuron: RS Response'));
xlabel('Time (in ms)');
ylabel('Spikes');

%RESULT: MANY SPIKES OBSERVED (68) DUE TO INCREASE IN MEAN SYNAPSE WEIGHT

end
function [syn_strength, num_iter] = problem_4(m_syn_stimulus,V_t_input, w0, sigma_w)
gamma = 1;
I0 = 1e-12;
tau = 15e-3;
tau_s = tau/4;
Ns = 100;
T = 0.5;
dt = 1e-4;
l_t = 1;
lambda = l_t*dt;
syn_strength = w0 + sigma_w*randn(Ns,1); %Strength gaussian -- randn -- standard normal
time = (dt:dt:T);
t_size = length(time);
num_iter = 0;
V_t = V_t_input; % to satisfy while condition
%[synap, time_of_stim ] = find(m_syn_stimulus == 1);
while(max(V_t) > -0.04)
    Iapplied = zeros(Ns,t_size);
for i = 1:t_size
    if(V_t(i) > -0.04) %Neuron 
        for t = 1:i-1
        j = find(m_syn_stimulus == t);
        if(isempty(j) == 1)
            continue;
        else
        del_spike = i-t;
        del_wk = -1*gamma*syn_strength(j,1)*(exp(-del_spike/tau) - exp(-del_spike/tau_s));
        syn_strength(j,1) = del_wk + syn_strength(j,1); 
        Iapplied(j,i) = Iapplied(j,i) + I0*syn_strength(j,1)*(exp(-(i - t)*dt/tau) - exp(-(i - t)*dt/tau_s));  
        plot_diag_wk(j, num_iter + 1) = del_wk;
        plot_diag_tk(j, num_iter + 1) = del_spike;
        break;
        end
        end
    end
end

% for k = 1: Ns
%     for i = 1:t_size
%         for t = 1:i
%             if (m_syn_stimulus(k) == t)
%                 Iapplied(k,i) = Iapplied(k,i) + I0*syn_strength(k,1)*(exp(-(i - t)*dt/tau) - exp(-(i - t)*dt/tau_s));  
%             end
%         end
%     end
% end 
I_total = sum(Iapplied, 1); % Sum of all synapses 1 x 5000
%Membrane Potential (AEF RS) :HW 1 Q3
V_t = euler_method(1,I_total,dt); %Same function as used in HW1
num_iter = num_iter + 1;
end
figure(6);
plot(time.*1000, V_t, 'linewidth', 2);
title(sprintf('Neuron: RS Response - No Spikes (4(a))'));
xlabel('Time (in ms)');
ylabel('Spikes');

%4(b)
figure(7);
for j = 1:num_iter
scatter(plot_diag_tk(:,j)*1000, plot_diag_wk(:,j));
hold on;
end
title('del_wk vs del_tk (4(b))');
xlabel('del_tk (in ms)');
ylabel('del_wk');
end
function problem_3(m_syn_stimulus)
gamma = 1;
w0 = 50;
sigma_w = 5;
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
V_t = ones(1, t_size)*-0.09; % to satisfy while condition
%plot_diag - Ns x 1 x 1 x num_iter
col = ['r', 'g', 'b', 'y', 'm', 'c', 'w', 'k']
while(max(V_t) < -0.04)
Iapplied = zeros(Ns,t_size);
for k = 1: Ns
    for i = 1:t_size
        for t = 1:i
            if (m_syn_stimulus(k,t) == 1)
                Iapplied(k,i) = Iapplied(k,i) + I0*syn_strength(k,1)*(exp(-(i - t)*dt/tau) - exp(-(i - t)*dt/tau_s));  
            end
        end
    end
end 
I_total = sum(Iapplied, 1); % Sum of all synapses 1 x 5000
%Membrane Potential (AEF RS) :HW 1 Q3
V_t = euler_method(1,I_total,dt); %Same function as used in HW1
[Max_value,Index] = max(V_t);
t_max = Index*dt;
%Update of Weights
for j = 1:Ns
    n = (Index - 1);
    flag = 0;
    while(n>0)
        if (m_syn_stimulus(j,n) == 1)
          flag = 1;
          flag_stimulus = n;
          break;
        end
        n = n-1;
    end
    if(flag == 1)
       del_t_k = (t_max - flag_stimulus*dt);
       del_wk = syn_strength(j,1)*(exp(-del_t_k/tau) - exp(-del_t_k/tau_s));
       syn_strength(j,1) = del_wk + syn_strength(j,1);
       plot_diag_wk(j, num_iter + 1) = del_wk;
       plot_diag_tk(j, num_iter + 1) = abs(del_t_k);
    end
end
num_iter = num_iter + 1;
end

fprintf('Number of Iterations: %d\n', num_iter);
fprintf('Weights: \n');
disp(syn_strength);

figure(4);
plot(time.*1000, V_t, 'linewidth', 2);
title(sprintf('Neuron: RS Response - First Single Spike (3(a))'));
xlabel('Time (in ms)');
ylabel('Spikes');

figure(5);
for j = 1:num_iter
    scatter(plot_diag_tk(:,j)*1000, plot_diag_wk(:,j));
    hold on;
end
title('del_wk vs del_tk (3(b))');
xlabel('del_tk (in ms)');
ylabel('del_wk');

end
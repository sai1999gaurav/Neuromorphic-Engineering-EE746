%% Problem 1
% PART (a)
T = 0.5;
dt = 0.1e-3;
k = 10;
lambda = k*dt;
time = (dt:dt:T);
poisson_output = makedist('Poisson','lambda',lambda);
t_size = length(time);
stimulus = zeros(1,t_size);
fprintf('Time Instants: ');
for i = 1:t_size
    spike = random(poisson_output,1);
    if (spike >= 1)
        stimulus(1,i) = 1; %spike-train
        fprintf('%f, ', i*dt); %timestamp
    end
end
figure(1);
subplot(3,1,1);
plot(time.*1000,stimulus(1,:), 'Linewidth',2);
title('Poisson Stimulus');
xlabel('Spikes');
ylabel('Time (in ms) ');
% PART (b) - Total Current Calculation
I0 = 1e-12;
we = 500;
tau = 15e-3;
tau_s = tau/4;
Iapplied = zeros(1,t_size);
for i = 1:t_size
    for t = 1:i
        if (stimulus(1,t) == 1)
            Iapplied(1,i) = Iapplied(1,i) + I0*we*(exp(-(i - t)*dt/tau) - exp(-(i - t)*dt/tau_s));  
        end
    end
end
subplot(3,1,2);
plot(time.*1000,Iapplied(1,:), 'Linewidth',2);
title('Synaptic Current');
ylabel('Iapplied');
xlabel('Time (in ms)');
%Membrane Potential (AEF RS) :HW 1 Q3
V_t = euler_method(1,Iapplied,dt); %Same function as used in HW1
subplot(3,1,3);
plot(time.*1000, V_t, 'linewidth', 2);
title(sprintf('Neuron: RS Response'));
xlabel('Time (in ms)');
ylabel('Spikes');
%% Problem 2a
Ns = 100;
w0 = 50;
sigma_w = 5;
T = 0.5;
dt = 1e-4;
l_t = 1;
lambda = l_t*dt;
syn_strength = w0 + sigma_w*randn(Ns,1); %Strength gaussian -- randn -- standard normal
time = (dt:dt:T);
poisson_output = makedist('Poisson', 'lambda', lambda);
t_size = length(time);
m_syn_stimulus = zeros(Ns, t_size);
for i = 1:Ns % Each Synapse
    for j = 1:t_size % Each time
        spike = random(poisson_output,1);
        if (spike >= 1)
         m_syn_stimulus(i,j) = 1; %spike-train
        end
    end
end
I0 = 1e-12;
tau = 15e-3;
tau_s = tau/4;
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
figure(2);
sgtitle(sprintf('Neuron Response with w0 = 50 and sigma_w = 5'));
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

%RESULT: NO SPIKES ISSUED, DUE TO LESS WEIGHTS

%PART (b) -- SAME STIMULUS
w0 = 250;
sigma_w = 25;
syn_strength = w0 + sigma_w*randn(Ns,1); %Strength gaussian -- randn -- standard normal
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
figure(3);
sgtitle(sprintf('Neuron Response with w0 = 250 and sigma_w = 25'));
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

% Problem 3 -- SAME STIMULUS
gamma = 1;
w0 = 50;
sigma_w = 5;
syn_strength = w0 + sigma_w*randn(Ns,1); %Strength gaussian -- randn -- standard normal
num_iter = 0;
V_t = ones(1, t_size)*-0.09; % to satisfy while condition
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
       syn_strength(j,1) = syn_strength(j,1)*(1 + exp(-del_t_k/tau) - exp(-del_t_k/tau_s));
    end
end
num_iter = num_iter + 1;
end

fprintf('Number of Iterations: %d\n', num_iter);
fprintf('Weights: ');
disp(syn_strength);

figure(4);
plot(time.*1000, V_t, 'linewidth', 2);
title(sprintf('Neuron: RS Response'));
xlabel('Time (in ms)');
ylabel('Spikes');

% %Result: 
% Number of Iterations: 5
% Weights:   165.2936
%    53.0303
%    47.3463
%    49.7780
%    49.3572
%    98.8588
%    45.7503
%    51.0064
%    52.7484
%    45.8921
%    47.0531
%    53.9118
%    40.2040
%    47.5685
%    51.6732
%    47.4775
%    48.7321
%    59.4139
%    42.5780
%    43.3906
%    49.7650
%    58.6208
%    48.4452
%    48.9178
%    53.7738
%    46.6784
%    54.5104
%    54.8192
%    52.5851
%    58.9157
%    42.2850
%    52.7052
%    49.8447
%    50.0257
%   123.2671
%    55.2158
%    44.1976
%    53.7434
%    81.7282
%    53.0447
%    53.6514
%    50.1975
%    45.2515
%   305.6106
%    51.9156
%    36.4064
%    50.7554
%    70.9463
%    45.9168
%    49.8506
%    48.9504
%    54.0080
%    49.6984
%    53.5643
%    48.9612
%    49.5255
%    40.8129
%    48.8495
%    51.4949
%   156.7559
%    56.4243
%    47.7855
%    41.5651
%    55.8315
%    58.2684
%    52.3795
%   185.7696
%    46.3691
%    58.2283
%    54.8906
%    53.8504
%    43.8901
%    56.5156
%    43.4044
%    46.6869
%    53.3584
%    50.6972
%    44.2660
%    44.8037
%    51.2461
%    50.6548
%    51.9624
%    71.0285
%    44.1940
%   108.0509
%    57.8358
%    50.4191
%    47.5861
%    51.6630
%    77.1421
%    45.3869
%    41.2184
%    46.1111
%    54.8149
%    46.7212
%    49.9354
%    45.0513
%    41.8302
%    46.3724
%    52.7744


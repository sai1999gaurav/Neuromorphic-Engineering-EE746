function problem_1()
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
sgtitle(sprintf('Problem 1'));
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
end
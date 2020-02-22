function problem_4(m_syn_stimulus)
gamma = 1;
w0 = 250;
sigma_w = 25;
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
for i = 1:t_size
    if(V_t(i) > -0.04) %Neuron Spiked
        
    end
end
end
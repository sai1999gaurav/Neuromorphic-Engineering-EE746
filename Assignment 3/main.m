global t1 t2 t3 T N M ;
Weight = [3000;3000;3000;3000;3000;3000];
Delay  = [1e-3;8e-3;5e-3;5e-3;9e-3;1e-3];

T = 100e-3;
del_t = 0.05e-3;
time = (del_t:del_t:T);
t_size = length(time); %size of time array
Iapp = 50e-9;
I_in = zeros(3,t_size);
%Case-1
for i = 1:9
    I_in(1,i) = Iapp;
    I_in(2,40+i) = Iapp;
    I_in(3,80+i) = Iapp;
end
V_t_pre = euler_method(I_in,del_t);
%hold on
%plot(time.*1000,V_t_pre(1,:), 'Linewidth',1);
%plot(time.*1000,V_t_pre(2,:), 'Linewidth',1);
%plot(time.*1000,V_t_pre(3,:), 'Linewidth',1);
%hold off

I_syn = zeros(6,t_size);
I_o = 1e-12;
tau = 15e-3/del_t;
tau_s = tau/4;
h = [t1;t1;t2;t2;t3;t3];
H =zeros(6,1);
Delay = Delay/del_t;
H(:,1) = h(:,1)+Delay(:,1);
for i = 1:t_size
    I_syn(1,i) =(I_o*Weight(1,1)*(exp(-(i-H(1,1))/tau)-exp(-(i-H(1,1))/tau_s)))*heaviside(i-H(1,1));
    I_syn(2,i) =(I_o*Weight(2,1)*(exp(-(i-H(2,1))/tau)-exp(-(i-H(2,1))/tau_s)))*heaviside(i-H(2,1));    
    I_syn(3,i) =(I_o*Weight(3,1)*(exp(-(i-H(3,1))/tau)-exp(-(i-H(3,1))/tau_s)))*heaviside(i-H(3,1));
    I_syn(4,i) =(I_o*Weight(4,1)*(exp(-(i-H(4,1))/tau)-exp(-(i-H(4,1))/tau_s)))*heaviside(i-H(4,1));
    I_syn(5,i) =(I_o*Weight(5,1)*(exp(-(i-H(5,1))/tau)-exp(-(i-H(5,1))/tau_s)))*heaviside(i-H(5,1));
    I_syn(6,i) =(I_o*Weight(6,1)*(exp(-(i-H(6,1))/tau)-exp(-(i-H(6,1))/tau_s)))*heaviside(i-H(6,1));
end
I_in_2 = zeros(3,t_size);
for i = 1:t_size
I_in_2(1,i) = I_syn(1,i)+I_syn(3,i)+I_syn(5,i);
I_in_2(2,i) = I_syn(2,i)+I_syn(4,i)+I_syn(6,i);
end
hold on
%plot(time.*1000,I_in(1,:), 'Linewidth',1);
%plot(time.*1000,I_in(2,:), 'Linewidth',1);
%plot(time.*1000,I_in(3,:), 'Linewidth',1);
%plot(time.*1000,I_syn(1,:), 'Linewidth',1);
%plot(time.*1000,I_syn(2,:), 'Linewidth',1);
%plot(time.*1000,I_syn(3,:), 'Linewidth',1);
%plot(time.*1000,I_syn(4,:), 'Linewidth',1);
%plot(time.*1000,I_syn(5,:), 'Linewidth',1);
%plot(time.*1000,I_syn(6,:), 'Linewidth',1);
%plot(time.*1000,heaviside(time.*1000-H(1,1)), 'Linewidth',1);
%plot(time.*1000,I_in_2(1,:), 'Linewidth',1);
hold off

V_t_post = euler_method(I_in_2,del_t);
hold on
subplot(2,1,1);
hold on
plot(time.*1000,I_syn(2,:), 'Linewidth',1);
plot(time.*1000,I_syn(4,:), 'Linewidth',1);
plot(time.*1000,I_syn(6,:), 'Linewidth',1);
plot(time.*1000,I_in_2(2,:), 'Linewidth',1);
title('Total Synaptic Current: Neuron b');
ylabel('I_{synapse}');
xlabel('Time (in ms)');
hold off
subplot(2,1,2);
plot(time.*1000,V_t_post(2,:), 'Linewidth',1);
title('Postsynaptic Potential: Neuron b');
ylabel('V_t');
xlabel('Time (in ms)');

%plot(time.*1000,V_t_post(2,:), 'Linewidth',1);
hold off
%%
Fanout = [1;5;1;5;1;5];
Weight = [3000;3000;3000;3000;3000;3000];
Delay  = [1e-3;8e-3;5e-3;5e-3;9e-3;1e-3];
T = 100e-3;
del_t = 0.05e-3;
time = (del_t:del_t:T);
t_size = length(time); %size of time array
Iapp = 50e-9;
I_in = zeros(3,t_size);
%Case-2
for i = 1:9
    I_in(3,i) = Iapp;
    I_in(2,30+i) = Iapp;
    I_in(1,70+i) = Iapp;
end
V_t_pre = euler_method(I_in,del_t);
%hold on
%plot(time.*1000,V_t_pre(1,:), 'Linewidth',1);
%plot(time.*1000,V_t_pre(2,:), 'Linewidth',1);
%plot(time.*1000,V_t_pre(3,:), 'Linewidth',1);
%hold off

I_syn = zeros(6,t_size);
I_o = 1e-12;
tau = 15e-3/del_t;
tau_s = tau/4;
h = [t1;t1;t2;t2;t3;t3];
H =zeros(6,1);
Delay = Delay/del_t;
H(:,1) = h(:,1)+Delay(:,1);
for i = 1:t_size
    I_syn(1,i) =(I_o*Weight(1,1)*(exp(-(i-H(1,1))/tau)-exp(-(i-H(1,1))/tau_s)))*heaviside(i-H(1,1));
    I_syn(2,i) =(I_o*Weight(2,1)*(exp(-(i-H(2,1))/tau)-exp(-(i-H(2,1))/tau_s)))*heaviside(i-H(2,1));    
    I_syn(3,i) =(I_o*Weight(3,1)*(exp(-(i-H(3,1))/tau)-exp(-(i-H(3,1))/tau_s)))*heaviside(i-H(3,1));
    I_syn(4,i) =(I_o*Weight(4,1)*(exp(-(i-H(4,1))/tau)-exp(-(i-H(4,1))/tau_s)))*heaviside(i-H(4,1));
    I_syn(5,i) =(I_o*Weight(5,1)*(exp(-(i-H(5,1))/tau)-exp(-(i-H(5,1))/tau_s)))*heaviside(i-H(5,1));
    I_syn(6,i) =(I_o*Weight(6,1)*(exp(-(i-H(6,1))/tau)-exp(-(i-H(6,1))/tau_s)))*heaviside(i-H(6,1));
end
I_in_2 = zeros(3,t_size);
for i = 1:t_size
I_in_2(1,i) = I_syn(1,i)+I_syn(3,i)+I_syn(5,i);
I_in_2(2,i) = I_syn(2,i)+I_syn(4,i)+I_syn(6,i);
end
%hold on
%plot(time.*1000,I_in(1,:), 'Linewidth',1);
%plot(time.*1000,I_in(2,:), 'Linewidth',1);
%plot(time.*1000,I_in(3,:), 'Linewidth',1);
%plot(time.*1000,I_syn(1,:), 'Linewidth',1);
%plot(time.*1000,I_syn(2,:), 'Linewidth',1);
%plot(time.*1000,I_syn(3,:), 'Linewidth',1);
%plot(time.*1000,I_syn(4,:), 'Linewidth',1);
%plot(time.*1000,I_syn(5,:), 'Linewidth',1);
%plot(time.*1000,I_syn(6,:), 'Linewidth',1);
%plot(time.*1000,heaviside(time.*1000-H(1,1)), 'Linewidth',1);
%plot(time.*1000,I_in_2(2,:), 'Linewidth',1);
%hold off

V_t_post = euler_method(I_in_2,del_t);
hold on
subplot(2,1,1);
hold on
plot(time.*1000,I_syn(1,:), 'Linewidth',1);
plot(time.*1000,I_syn(3,:), 'Linewidth',1);
plot(time.*1000,I_syn(5,:), 'Linewidth',1);
plot(time.*1000,I_in_2(1,:), 'Linewidth',1);
title('Total Synaptic Current: Neuron a');
ylabel('I_{synapse}');
xlabel('Time (in ms)');
hold off
subplot(2,1,2);
plot(time.*1000,V_t_post(1,:), 'Linewidth',1);
title('Postsynaptic Potential: Neuron a');
ylabel('V_t');
xlabel('Time (in ms)');
hold off
%%
%Q2
N = 500; %total number of neurons
M = 25; %excited neurons
%first 400 neurons are excitatory and last 100 are inhibitory
for i=1:N
    spike_time{i} = [];
    arrival_time{i} = [];
    strength{i} = [];
    pre_neuron{i} = [];
end

T = 1000e-3;
del_t = 1e-3;
time = (del_t:del_t:T);
t_size = length(time); %size of time array
l_t = 100;
lambda = l_t*del_t;
poisson_output = makedist('Poisson', 'lambda', lambda);
m_syn_stimulus = zeros(M, t_size);
for i = 1:M % Each Synapse
    for j = 1:t_size % Each time
        spike = random(poisson_output,1);
        if (spike >= 1)
         m_syn_stimulus(i,j) = 1; %spike-train
        end
    end
end
I_o = 1e-12;
w_e = 3000;
tau = (1+19*rand)*1e-3;
tau_s = tau/4;
I = zeros(N,t_size);
for k = 1: M
    for i = 1:t_size
        for t = 1:i
            if (m_syn_stimulus(k,t) == 1)
                I(k,i) = I(k,i) + I_o*w_e*(exp(-(i - t)*del_t/tau) - exp(-(i - t)*del_t/tau_s));  
            end
        end
    end
end 
%{
figure(2);
%sgtitle(sprintf('Neuron Response with w0 = 50 and sigma_w = 5 (2(a))'));
subplot(2,1,1);
plot(time.*1000,Iapplied(1,:), 'Linewidth',2);
title('Input Current');
ylabel('Iapplied');
xlabel('Time (in ms)');
%}
V_t = zeros(N,t_size);
[V_t,spike_time] = euler_method_general(I,del_t,spike_time);
subplot(2,1,2);
plot(time.*1000, V_t(1,:), 'linewidth', 2);
title(sprintf('Neuron Response(V)'));
xlabel('Time (in ms)');
ylabel('Spikes');


%giving synaptic current with delay
for i=1:N
    Fanout{i}=[];
end    
for i = 1:(N-100)
    r = randi([1 N],1,N/10);
    Fanout{i} = [Fanout{i} r];
end    
for i = (N-99):N
    s = randi([1 N-100],1,N/10);
    Fanout{i} = [Fanout{i} s];
end   
%updating arrival time, strength, pre_neuron
for i = 1:M    %number of initial neurons given input
    for k=1:length(spike_time{i})   %total number of spikes in an iniital neuron
            for j=1:N/10    %FanOut of every initial neuron
                temp = 1;   %inhibitory delay
                temp1 = randi([1,20],1,1);  %exhibitory delay
                a = Fanout{i}(1,j); 
                if(a<=400) %excitatory
                    temp2 = spike_time{i}(1,k);
                    temp3 = {temp1+temp2};
                    arrival_time{a} =  [arrival_time{a} temp3];
                    strength{a} = [strength{a} w_e];
                    pre_neuron{a} = [pre_neuron{a} i];
                end
                    temp2 = spike_time{i}(1,k);
                    temp3 = {temp+temp2};
                    arrival_time{a} =  [arrival_time{a} temp3]; %inhibitory
                    strength{a} = [strength{a} -w_e];
                    pre_neuron{a} = [pre_neuron{a} i];
            end
    end
end
temporary = cell([N,3]);
 for i=1:N
    temporary{i,1} = (transpose(arrival_time{i}));
    temporary{i,2} = (transpose(num2cell(strength{i})));
    temporary{i,3} = (transpose(num2cell(pre_neuron{i})));
 end
temporary2 = cell([N,1]);
temporary3 = cell([N,1]);
 for i=1:N
     temporary2{i,1} = cell2mat([temporary{i,1} temporary{i,2} temporary{i,3}]);
     temporary3{i,1} = sortrows(temporary2{i,1},1); 
 end    
for i=1:N %current in excitatory neurons
    for j=1:t_size
      if(temporary3{i,1}())  
     I(i,j) = I(i,j) + I_o*w_e*(exp(-(i - t)*del_t/tau) - exp(-(i - t)*del_t/tau_s));
    end
end    

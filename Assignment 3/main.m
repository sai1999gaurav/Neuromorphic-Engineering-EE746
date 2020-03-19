global t1 t2 t3 ;
Fanout = [1;5;1;5;1;5];
Weight = [3000;3000;3000;3000;3000;3000];
Delay  = [1e-3;8e-3;5e-3;5e-3;9e-3;1e-3];

for i=1:5
    spike_time{i} = [];
    arrival_time{i} = [];
    strength{i} = [];
    pre_neuron{i} = [];
end

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

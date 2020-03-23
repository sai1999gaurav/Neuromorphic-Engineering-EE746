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
    Fanout{i}=[];    
    spike_time_final{i} = [];
    arrival_time{i} = [];
    strength{i} = [];
    pre_neuron{i} = [];
end

T = 1000e-3;
del_t = 0.25e-3;
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
tau = (randi([1,20]))*1e-3/del_t;
tau_s = tau/4;
I = zeros(N,t_size);
for k = 1: M
    for i = 1:t_size
        for t = 1:i
            if (m_syn_stimulus(k,t) == 1)
                I(k,i) = I(k,i) + I_o*w_e*(exp(-(i - t)/tau) - exp(-(i - t)/tau_s));  
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
M=500;
[V_t,spike_time] = euler_method_general(I,del_t,spike_time,M);
M=25;
subplot(2,1,2);
plot(time.*1000, V_t(1,:), 'linewidth', 2);
title(sprintf('Neuron Response(V)'));
xlabel('Time (in ms)');
ylabel('Spikes');


%giving synaptic current with delay   
for i = 1:(N-100)
    r = randi([1 N],1,N/10);
    Fanout{1,i} = [Fanout{1,i} r];
end    
for i = (N-99):N
    s = randi([1 N-100],1,N/10);
    Fanout{1,i} = [Fanout{1,i} s];
end   

%updating arrival time, strength, pre_neuron
for i = 1:N    %number of initial neurons given input
    null_spike = isempty(spike_time{i});
    if (null_spike ==1)
        continue;
    end    
    for k=1:length(spike_time{i})   %total number of spikes in an initial neuron
            for j=1:N/10    %FanOut of every initial neuron
                temp = 1;   %inhibitory delay
                temp1 = randi([1,20],1,1);  %exhibitory delay
                a = Fanout{i}(1,j); 
                if(i<=400) %excitatory
                    temp2 = spike_time{i}(1,k);
                    temp3 = {temp1+temp2};
                    arrival_time{a} =  [arrival_time{a} temp3];
                    strength{a} = [strength{a} w_e];
                    pre_neuron{a} = [pre_neuron{a} i];
                end
                if(i>400)   %inhibitory
                    temp2 = spike_time{i}(1,k);
                    temp3 = {temp+temp2};
                    arrival_time{a} =  [arrival_time{a} temp3]; %inhibitory
                    strength{a} = [strength{a} -w_e];
                    pre_neuron{a} = [pre_neuron{a} i];
                end   
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
     null_matrix = isempty(temporary2{i,1});
     if(null_matrix == 1) %skip the iteration if it's a null matrix
        continue;
     end
     temporary3{i,1} = sortrows(temporary2{i,1},1); 
 end
 t_size = T/del_t;
 I_2 = zeros(N,t_size); %Second iteration of current    
 I_o = 1e-12;
 tau = (randi([1,20]))*1e-3/del_t;
 tau_s = tau/4;
 for i=1:N %current in excitatory neurons
      [m,n] = size(temporary3{i,1});    %to get size of arrival matrix, n will be 3 always
      for j=1:m 
            arrivaltime = temporary3{i,1}(j,1);
            weight = temporary3{i,1}(j,2);
            null_flag = isempty(arrivaltime);
              if(arrivaltime>T*1e3) || (null_flag == 1)  %skip iteration if arrival time is more than t_size
                  continue;
              end
            for k=1:t_size
                if(k >= arrivaltime*1e-3/del_t)
                I_2(i,k) =I_2(i,k)+(I_o*(weight)*(exp(-(k - arrivaltime*1e-3/del_t)/tau) - exp(-(k - arrivaltime*1e-3/del_t)/tau_s)));
                end
            end
      end
 end
 %Total current 
 I_total = zeros(N,t_size);
 for i=1:N
    I_total(i,:) = I(i,:) + I_2(i,:);
 end
 M=500;
 [V_t,spike_time_final] = euler_method_general(I_total,del_t,spike_time_final,M);
 
 %{
 %second iteration
 for i = 1:N    %number of initial neurons given input
    null_spike = isempty(spike_time_final{i});
    if (null_spike == 1)
        continue;
    end    
    for k=1:length(spike_time_final{i})   %total number of spikes in an iniital neuron
            for j=1:N/10    %FanOut of every initial neuron
                temp = 1;   %inhibitory delay
                temp1 = randi([1,20],1,1);  %exhibitory delay
                a = Fanout{i}(1,j); 
                if(i<=400) %excitatory
                    temp2 = spike_time_final{i}(1,k);
                    temp3 = {temp1+temp2};
                    arrival_time{a} =  [arrival_time{a} temp3];
                    strength{a} = [strength{a} w_e];
                    pre_neuron{a} = [pre_neuron{a} i];
                end
                if(i>400)   %inhibitory
                    temp2 = spike_time_final{i}(1,k);
                    temp3 = {temp+temp2};
                    arrival_time{a} =  [arrival_time{a} temp3]; %inhibitory
                    strength{a} = [strength{a} -w_e];
                    pre_neuron{a} = [pre_neuron{a} i];
                end   
            end
    end
end
%temporary = cell([N,3]);
 for i=1:N
    temporary{i,1} = (transpose(arrival_time{i}));
    temporary{i,2} = (transpose(num2cell(strength{i})));
    temporary{i,3} = (transpose(num2cell(pre_neuron{i})));
 end
%temporary2 = cell([N,1]);
%temporary3 = cell([N,1]);
 for i=1:N
     temporary2{i,1} = cell2mat([temporary{i,1} temporary{i,2} temporary{i,3}]);
     null_matrix = isempty(temporary2{i,1});
     if(null_matrix == 1) %skip the iteration if it's a null matrix
        continue;
     end
     temporary3{i,1} = sortrows(temporary2{i,1},1); 
 end
 t_size = T/del_t;
 I_3 = zeros(N,t_size); %Second iteration of current    
 I_o = 1e-12;
 tau = (randi([1,20]))*1e-3/del_t;
 tau_s = tau/4;
 for i=1:N %current in excitatory neurons
      [m,n] = size(temporary3{i,1});    %to get size of arrival matrix, n will be 3 always
      for j=1:m 
            arrivaltime = temporary3{i,1}(j,1);
            weight = temporary3{i,1}(j,2);
            null_flag = isempty(arrivaltime);
              if(arrivaltime>T*1e3) || (null_flag == 1)  %skip iteration if arrival time is more than t_size
                  continue;
              end
            for k=1:t_size
                if(k >= arrivaltime*1e-3/del_t)
                I_3(i,k) = I_3(i,k)+(I_o*(weight)*(exp(-(k - arrivaltime*1e-3/del_t)/tau) - exp(-(k - arrivaltime*1e-3/del_t)/tau_s)));
                end
            end
      end
 end
 %Total current 
 I_total = zeros(N,t_size);
 for i=1:N
    I_total(i,:) = I(i,:) + I_2(i,:) + I_3(i,:);
 end
 M=500;
 [V_t,spike_time_final] = euler_method_general(I_total,del_t,spike_time_final,M);
 %}
 %rastor plot, spiking information is in spike_time cell array
 RastorPlot = zeros(N,t_size);
 for i = 1:N
    for j = 1:length(spike_time_final{1,i})
        RastorPlot(i,int32(spike_time_final{1,i}(1,j)*1e-3/del_t)) = 1;
    end    
 end
%{
for j = 1:20
hold on
figure(j+2)
subplot(25,1,1);
plot(time.*1000,RastorPlot((1+25*(j-1)),:), 'Linewidth',0.5);
set(gca,'XTick',[],'YTick',[]);
title(['Neuron' num2str((1+25*(j-1))) '-' num2str(25*j)]);
for i = 2:24    
        subplot(25,1,i);
        plot(time.*1000,RastorPlot(25*(j-1)+i,:), 'Linewidth',0.5);
        set(gca,'XTick',[],'YTick',[]);
end
subplot(25,1,25);
plot(time.*1000,RastorPlot(25*j,:), 'Linewidth',0.5);
set(gca,'Ytick',[]);
xlabel('Time(in ms)');
hold off
end
%}
%Q2b

figure(20)
hold on
c = linspace(1,N,N);
grid on
grid minor 
for k = 1:t_size    
    for i = 1:N
    if(RastorPlot(i,k) == 0)
        continue;
    end    
        scatter(k,RastorPlot(i,k)+i-1,1,c(1,i),'filled');
    end    
end    
xlabel('Time(in ms)');
ylabel('Spikes');
hold off

%printing fanout
figure(21);
hold on
c = linspace(1,500,500);
grid on
grid minor 
for i = 1:N
    for k = 1:N/10    
        scatter(i,Fanout{1,i}(1,k),1,c(1,i),'filled');
    end    
end    
xlabel('Neurons');
ylabel('FanOut');
hold off

spike_time_all = zeros(2,t_size);
for i = 1:N
    if(i<=400)
        null_spike_matrix_flag = isempty(spike_time_final{1,i});
        if (null_spike_matrix_flag == 1)
            continue;
        end    
        for j = 1:length(spike_time_final{1,i})
            spike_time_all(1,int32(spike_time_final{1,i}(1,j)*1e-3/del_t)) = 1;
        end
    end    
    if (i>400)   
        null_spike_matrix_flag = isempty(spike_time_final{1,i});
        if (null_spike_matrix_flag == 1)
            continue;
        end
        for j = 1:length(spike_time_final{1,i})
            spike_time_all(2,int32(spike_time_final{1,i}(1,j)*1e-3/del_t)) = 1;
        end        
    end   
end    
R_e = zeros(1,t_size);
R_i = zeros(1,t_size);
window = 10e-3/del_t;

for i=1:t_size-window
    for j=1:window+1
        R_e(1,i) = R_e(1,i)+spike_time_all(1,i+j-1);
        R_i(1,i) = R_i(1,i)+spike_time_all(2,i+j-1);
    end    
end    

figure(22)
hold on
    plot(time/del_t,R_e(1,:));
    plot(time/del_t,R_i(1,:));
    legend('excitatory neurons','inhibitory neurons');
hold off




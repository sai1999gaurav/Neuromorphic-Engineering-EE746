global T N M ;
del_t = 0.05e-3;
Fanout = {[]; [1 5]; [1 5]; [1 5]; []};
Weight = {[]; [3000 3000]; [3000 3000]; [3000 3000]; []};
Delay  = {[]; [1e-3/del_t 8e-3/del_t]; [5e-3/del_t 5e-3/del_t]; [9e-3/del_t 1e-3/del_t]; []};

T = 100e-3;
time = (del_t:del_t:T);
t_size = length(time); %size of time array
Iapp = 50e-9;
I_in = zeros(5,t_size);
pulse_width = 1e-3/del_t;
for i=1:5
    spike_time{i,1} = [];    
    spike_time_final{i,1} = [];
    arrival_time{i,1} = [];
    pre_neuron{i,1} = [];
end

%Case-1
for i = 1:pulse_width
    I_in(4,i) = Iapp;
    I_in(3,(3e-3/del_t)+i) = Iapp;
    I_in(2,(7e-3/del_t)+i) = Iapp;
end
[V_t_pre,spike_time] = euler_method(I_in,del_t,spike_time);
%updating arrival time and pre-neuron
for i=1:5
    Fanout_empty_flag = isempty(Fanout{i,1});
    if (Fanout_empty_flag == 1)
        continue;
    end    
    for j=1:length(spike_time{i,1})
        for k = 1:length(Fanout{i,1})
            arrival_time_var = spike_time{i,1}(1,j)+Delay{i,1}(1,k);
            arrival_time{Fanout{i,1}(1,k),1} = [arrival_time{Fanout{i,1}(1,k),1} arrival_time_var];
            pre_neuron{Fanout{i,1}(1,k),1} = [pre_neuron{Fanout{i,1}(1,k),1} i];
        end
    end
end    
figure(1)
hold on
subplot(2,1,1);
plot(time.*1000,I_in(2,:),'DisplayName','neuron-b', 'Linewidth',1);
title('Current given to neuron-b');
xlabel('Time(in ms)');
ylabel('I');
subplot(2,1,2);
plot(time.*1000,V_t_pre(2,:),'DisplayName','neuron-b', 'Linewidth',1);
title('Potential generated in neuron-b');
xlabel('Time(in ms)');
ylabel('V');
hold off
figure(2)
hold on
subplot(2,1,1);
plot(time.*1000,I_in(3,:),'DisplayName','neuron-c', 'Linewidth',1);
title('Current given to neuron-c');
xlabel('Time(in ms)');
ylabel('I');
subplot(2,1,2);
plot(time.*1000,V_t_pre(3,:),'DisplayName','neuron-c', 'Linewidth',1);
title('Potential generated in neuron-c');
xlabel('Time(in ms)');
ylabel('V');
hold off
figure(3)
hold on
subplot(2,1,1);
plot(time.*1000,I_in(4,:),'DisplayName','neuron-d', 'Linewidth',1);
title('Current given to neuron-d');
xlabel('Time(in ms)');
ylabel('I');
subplot(2,1,2);
plot(time.*1000,V_t_pre(4,:),'DisplayName','neuron-d', 'Linewidth',1);
title('Potential generated in neuron-d');
xlabel('Time(in ms)');
ylabel('V');
hold off
I_syn = zeros(6,t_size);
I_o = 1e-12;
tau = 15e-3/del_t;
tau_s = tau/4;
w_e = 3000;
for j = [1 5]
        for i = 1:t_size
            if (j == 1)
                I_syn(1,i) = I_o*w_e*(exp(-(i-arrival_time{j,1}(1,1))/tau)-exp(-(i-arrival_time{j,1}(1,1))/tau_s))*heaviside(i-arrival_time{j,1}(1,1));
                I_syn(2,i) = I_o*w_e*(exp(-(i-arrival_time{j,1}(1,2))/tau)-exp(-(i-arrival_time{j,1}(1,2))/tau_s))*heaviside(i-arrival_time{j,1}(1,2));
                I_syn(3,i) = I_o*w_e*(exp(-(i-arrival_time{j,1}(1,3))/tau)-exp(-(i-arrival_time{j,1}(1,3))/tau_s))*heaviside(i-arrival_time{j,1}(1,3));
            end
            if (j == 5)
                I_syn(4,i) = I_o*w_e*(exp(-(i-arrival_time{j,1}(1,1))/tau)-exp(-(i-arrival_time{j,1}(1,1))/tau_s))*heaviside(i-arrival_time{j,1}(1,1));
                I_syn(5,i) = I_o*w_e*(exp(-(i-arrival_time{j,1}(1,2))/tau)-exp(-(i-arrival_time{j,1}(1,2))/tau_s))*heaviside(i-arrival_time{j,1}(1,2));
                I_syn(6,i) = I_o*w_e*(exp(-(i-arrival_time{j,1}(1,3))/tau)-exp(-(i-arrival_time{j,1}(1,3))/tau_s))*heaviside(i-arrival_time{j,1}(1,3));
            end
        end
end
figure(4)
hold on
%plot(time.*1000,V_t_pre(2,:),'DisplayName','neuron-b', 'Linewidth',1);
%plot(time.*1000,V_t_pre(3,:),'DisplayName','neuron-b', 'Linewidth',1);
%plot(time.*1000,V_t_pre(4,:),'DisplayName','neuron-b', 'Linewidth',1);
plot(time.*1000,I_syn(1,:), 'Linewidth',1);
plot(time.*1000,I_syn(2,:), 'Linewidth',1);
plot(time.*1000,I_syn(3,:), 'Linewidth',1);
plot(time.*1000,I_syn(4,:), 'Linewidth',1);
plot(time.*1000,I_syn(5,:), 'Linewidth',1);
plot(time.*1000,I_syn(6,:), 'Linewidth',1);
legend('synaptic-current-1(neuro:b-a)','synapse-current-2(neuro:c-a)','synapse-current-3(neuro:d-a)','synapse-current-4(neuro:b-e)','synapse-current-5(neuro:c-e)','synapse-current-6(neuro:d-e)');
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
I_post = zeros(5,t_size);
I_post(1,:)=I_syn(1,:)+I_syn(2,:)+I_syn(3,:);
I_post(5,:)=I_syn(4,:)+I_syn(5,:)+I_syn(6,:);
[V_t_post,spike_time_final] = euler_method(I_post,del_t,spike_time_final);
figure(5)
subplot(2,1,1);
hold on
plot(time.*1000,I_post(1,:), 'Linewidth',1);
plot(time.*1000,I_post(5,:), 'Linewidth',1);
legend('neuron-a','neuron-e');
hold off
title('Total Synaptic Current: Neuron a,e');
ylabel('I_{synapse}');
xlabel('Time (in ms)');
subplot(2,1,2);
hold on
plot(time.*1000,V_t_post(1,:), 'Linewidth',1);
plot(time.*1000,V_t_post(5,:), 'Linewidth',1);
legend('neuron-a','neuron-e');
hold off
title('Postsynaptic Potential: Neuron a,e');
ylabel('V_t');
xlabel('Time (in ms)');
%plot(time.*1000,V_t_post(2,:), 'Linewidth',1);

%%
%Q2
N = 500; %total number of neurons
E = (20/100)*N;
M = 25; %excited neurons
%first 400 neurons are excitatory and last 100 are inhibitory
for i=1:N
    spike_time{i} = [];    
    spike_time_2{i} = [];
    spike_time_3{i} = [];
    spike_time_4{i} = [];
    arrival_time{i} = [];
    arrival_time_2{i} = [];
    arrival_time_3{i} = [];
    strength{i} = [];
    strength_2{i} = [];
    strength_3{i} = [];
    pre_neuron{i} = [];
    pre_neuron_2{i} = [];
    pre_neuron_3{i} = [];
end
del_t = 0.5e-3;
%defining axonal delays
axonal_delay = zeros(1,N);
for i = 1:N-E
    axonal_delay(1,i) = randi([1,20],1,1)*1e-3/del_t; 
end    
for i = N-E+1:N
    axonal_delay(1,i) = 1e-3/del_t;
end    
T = 1000e-3;
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
tau = 15e-3/del_t;
tau_s = tau/4;
I = zeros(N,t_size);
%giving input Poisson stimulus to first 25 neurons
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
%finding spikes in first 25 neurons, first iteration
[V_t,spike_time] = euler_method_general(I,del_t,spike_time,N);
M=25;
figure(1)
subplot(2,1,2);
plot(time.*1000, V_t(1,:), 'linewidth', 2);
title(sprintf('Neuron Response(V)'));
xlabel('Time (in ms)');
ylabel('Spikes');
%%
%defining FanOut   
for i = 1:(N-E)
    r = randi([1 N],1,N/10+3);
    r = unique(r);
    Fanout{1,i} = [Fanout{1,i} r];
end    
for i = (N-E+1):N
    s = randi([1 N-E],1,N/10+3);
    s = unique(s);
    Fanout{1,i} = [Fanout{1,i} s];
end   
%%
%updating arrival time, strength, pre_neuron to all neurons
for i = 1:N    %number of initial neurons given input
    null_spike = isempty(spike_time{1,i}); %skip the iteration if a neuron doesn't get spikes
    if (null_spike ==1)
        continue;
    end    
    for k=1:length(spike_time{1,i})   %total number of spikes in an initially excited neuron
            for j=1:N/10    %FanOut of every initially excited neuron
                a = Fanout{1,i}(1,j); 
                if(i<=400) %excitatory
                    temp2 = spike_time{1,i}(1,k);
                    temp3 = axonal_delay(1,i)+temp2;
                    if (temp3*1e-3/del_t > t_size)
                        continue;
                    end    
                    arrival_time{a} =  [arrival_time{a} temp3];
                    strength{a} = [strength{a} w_e];
                    pre_neuron{a} = [pre_neuron{a} i];
                end
                if(i>400)   %inhibitory
                    temp2 = spike_time{1,i}(1,k);
                    temp3 = axonal_delay(1,i)+temp2;
                    if (temp3*1e-3/del_t > t_size)
                        continue;
                    end
                    arrival_time{a} =  [arrival_time{a} temp3]; %inhibitory
                    strength{a} = [strength{a} -w_e];
                    pre_neuron{a} = [pre_neuron{a} i];
                end   
            end
    end
end
%{
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
 %}
 t_size = T/del_t;
 I_2 = zeros(N,t_size);    
 I_o = 1e-12;
 tau = 15e-3/del_t;
 tau_s = tau/4;
 %finding current in all neurons after they receive their first spikes
 for i=1:N 
      [m,n] = size(arrival_time{1,i});    %to get size of arrival matrix, n will be 1 always
      if (m == 0)   %skip iteration if the neuron doesn't get any spikes
          continue;
      end    
      for j=1:n 
            arrivaltime = arrival_time{1,i}(1,j);
            weight = strength{1,i}(1,j);
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
 %finding spikes in all neurons, second iteration
 [V_t_2,spike_time_2] = euler_method_general(I_total,del_t,spike_time_2,M);
 %%
 %finding spikes in all excitatory and inhibitory neurons second time
 for i = 1:N    
    null_spike = isempty(spike_time_2{1,i});
    if (null_spike == 1)
        continue;
    end    
    for k=1:length(spike_time_2{1,i})   %total number of spikes in an iniital neuron
            for j=1:N/10    %FanOut of every initial neuron
                a = Fanout{i}(1,j); 
                if(i<=400) %excitatory
                    temp2 = spike_time_2{1,i}(1,k);
                    temp3 = axonal_delay(1,i)+temp2;
                    if (temp3*1e-3/del_t > t_size)
                        continue;
                    end
                    arrival_time_2{a} =  [arrival_time_2{a} temp3];
                    strength_2{a} = [strength_2{a} w_e];
                    pre_neuron_2{a} = [pre_neuron_2{a} i];
                end
                if(i>400)   %inhibitory
                    temp2 = spike_time_2{1,i}(1,k);
                    temp3 = axonal_delay(1,i)+temp2;
                    if (temp3*1e-3/del_t > t_size)
                        continue;
                    end
                    arrival_time_2{a} =  [arrival_time_2{a} temp3]; %inhibitory
                    strength_2{a} = [strength_2{a} -w_e];
                    pre_neuron_2{a} = [pre_neuron_2{a} i];
                end   
            end
    end
 end
%{
temporary_2 = cell([N,3]);
 for i=1:N
    temporary_2{i,1} = (transpose(arrival_time_2{i}));
    temporary_2{i,2} = (transpose(num2cell(strength_2{i})));
    temporary_2{i,3} = (transpose(num2cell(pre_neuron_2{i})));
 end
temporary2_2 = cell([N,1]);
temporary3_2 = cell([N,1]);
 for i=1:N
     temporary2_2{i,1} = cell2mat([temporary_2{i,1} temporary_2{i,2} temporary_2{i,3}]);
     null_matrix = isempty(temporary2_2{i,1});
     if(null_matrix == 1) %skip the iteration if it's a null matrix
        continue;
     end
     temporary3_2{i,1} = sortrows(temporary2_2{i,1},1); 
 end
 %}
 t_size = T/del_t;
 I_3 = zeros(N,t_size); %third iteration of current    
 I_o = 1e-12;
 tau = 15e-3/del_t;
 tau_s = tau/4;
 for i=1:N %current in excitatory neurons
      [m,n] = size(arrival_time_2{1,i});    %to get size of arrival matrix, n will be 1 always
      if (m ==0)    %if a neuron doesn't get any spikes then skip the iteration
          continue;
      end    
      for j=1:n 
            arrivaltime = arrival_time_2{1,i}(1,j);
            weight = strength_2{1,i}(1,j);
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
 for i=1:N
    I_total(i,:) = I_total(i,:) + I_3(i,:);
 end
 
 %finding spikes in all neurons, third iteration
 [V_t_3,spike_time_3] = euler_method_general(I_total,del_t,spike_time_3,N);
%%
 %fourth iteration
 for i = 1:N    
    null_spike = isempty(spike_time_3{1,i});
    if (null_spike == 1)
        continue;
    end    
    for k=1:length(spike_time_3{1,i})   
            for j=1:N/10    %FanOut of every neuron
                a = Fanout{i}(1,j); 
                if(i<=400) %excitatory
                    temp2 = spike_time_3{1,i}(1,k);
                    temp3 = axonal_delay(1,i)+temp2;
                    if (temp3*1e-3/del_t > t_size)
                        continue;
                    end
                    arrival_time_3{a} =  [arrival_time_3{a} temp3];
                    strength_3{a} = [strength_3{a} w_e];
                    pre_neuron_3{a} = [pre_neuron_3{a} i];
                end
                if(i>400)   %inhibitory
                    temp2 = spike_time_3{1,i}(1,k);
                    temp3 = axonal_delay(1,i)+temp2;
                    if (temp3*1e-3/del_t > t_size)
                        continue;
                    end
                    arrival_time_3{a} =  [arrival_time_3{a} temp3]; %inhibitory
                    strength_3{a} = [strength_3{a} -w_e];
                    pre_neuron_3{a} = [pre_neuron_3{a} i];
                end   
            end
    end
 end

 t_size = T/del_t;
 I_4 = zeros(N,t_size); %fourth iteration of current    
 I_o = 1e-12;
 tau = 15e-3/del_t;
 tau_s = tau/4;
 for i=1:N %current in excitatory neurons
      [m,n] = size(arrival_time_3{1,i});    %to get size of arrival matrix, n will be 1 always
      if (m ==0)    %if a neuron doesn't get any spikes then skip the iteration
          continue;
      end    
      for j=1:n 
            arrivaltime = arrival_time_3{1,i}(1,j);
            weight = strength_3{1,i}(1,j);
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
 for i=1:N
    I_total(i,:) = I_total(i,:) + I_4(i,:);
 end
 
 %finding spikes in all neurons, fourth iteration
 [V_t_4,spike_time_4] = euler_method_general(I_total,del_t,spike_time_4,N);
%%

 %rastor plot, spiking information is in spike_time cell array
 RastorPlot = zeros(N,t_size);
 for i = 1:N
    for j = 1:length(spike_time_4{1,i})
        RastorPlot(i,int32(spike_time_4{1,i}(1,j)*1e-3/del_t)) = 1;
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
%%
 
%Q2b
figure(2)
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

 %%

%printing fanout
figure(3);
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

 %%

spike_time_all = zeros(2,t_size);
for i = 1:N
    if(i<=400)
        null_spike_matrix_flag = isempty(spike_time_4{1,i});
        if (null_spike_matrix_flag == 1)
            continue;
        end    
        for j = 1:length(spike_time_4{1,i})
            spike_time_all(1,int32(spike_time_4{1,i}(1,j)*1e-3/del_t)) = 1;
        end
    end    
    if (i>400)   
        null_spike_matrix_flag = isempty(spike_time_4{1,i});
        if (null_spike_matrix_flag == 1)
            continue;
        end
        for j = 1:length(spike_time_4{1,i})
            spike_time_all(2,int32(spike_time_4{1,i}(1,j)*1e-3/del_t)) = 1;
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

figure(4)
hold on
    plot(time/del_t,R_e(1,:));
    plot(time/del_t,R_i(1,:));
    legend('excitatory neurons','inhibitory neurons');
hold off




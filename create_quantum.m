close all
clear all
clc

%% Stimulation signal script

% Create stimulation protocol_1 with a random order 
amp = [0.1, 0.3; 0.1, 0.45; 0.1, 0.6; 0.05, 0.6; 0.15, 0.2; 0.2, 0.6];
for x = 1:length(amp)
random_amp_pairs_1= amp(randperm(size(amp, 1)), :);
end

one_quantum = random_amp_pairs_1;

save('1_quantum', 'one_quantum');

%% Create waveforms 
% 1 electrode pair - quantum

I1x = [];

dt = 0.004;
carrier_f = 2000;
pulse_f = 10;

f1 = carrier_f/1000;
f2 = (carrier_f + pulse_f)/1000;

ramp_up_t = 500; % ms
ramp_down_t = 500;

for a = 1:length(random_amp_pairs_1)
        
A1 = random_amp_pairs_1(a,1);
A2 = random_amp_pairs_1(a,2);

phi1 = 0;
phi2 = pi;
    
each_stim_t = 6*1000; 
each_break_t = 5*1000;
total_t = 71*1000; 

stim_tt = dt:dt:each_stim_t;

I1_stim(a,:) = A1*cos(2*pi*f1*stim_tt+phi1) + A2*cos(2*pi*f2*stim_tt+phi1);

if ramp_up_t
    ramp_up_size = round(ramp_up_t/dt);
    ramp_vec = 0:1/ramp_up_size:1;

    I1_stim(a,1:length(ramp_vec)) = ramp_vec.*I1_stim(a,1:length(ramp_vec));

end

if ramp_down_t
    ramp_down_size = round(ramp_down_t/dt);
    ramp_vec = fliplr(0:1/ramp_down_size:1);

    I1_stim(a,end-length(ramp_vec)+1:end) = ramp_vec.*I1_stim(a,end-length(ramp_vec)+1:end);
    
end

% figure,plot(I1_stim(a,:))

% period of zeros

break_tt = dt:dt:each_break_t;

I1_break= zeros(1,length(break_tt));

% complete 'cycle'

I1_cycle(a,:) = [I1_stim(a,:), I1_break];

I1x = [I1x, I1_cycle(a,:)];
end

each_pre_t = 5*1000;
pre_tt = dt:dt:each_pre_t;
I1_pre = zeros(1,length(pre_tt));

I1x = [I1_pre, I1x];
I2x = zeros(size(I1x));
figure,plot(I1x + I2x)
% saveas(gcf,'quantum_1ep', '.fig');

I1 = I1x;
I2 = I2x;
save('quantum_1ep','I1','I2');

% %% Create waveforms 
% % 2 electrode pair - quantum
% 
% I1x = [];
% I2x = [];
% 
% dt = 0.004;
% carrier_f = 2000;
% pulse_f = 10;
% 
% f1 = carrier_f/1000;
% f2 = (carrier_f + pulse_f)/1000;
% 
% ramp_up_t = 500; % ms
% ramp_down_t = 500;
% 
% for a = 1:length(random_amp_pairs_1)
%         
% A1 = random_amp_pairs_1(a,1);
% A2 = random_amp_pairs_1(a,2);
% 
% phi1 = 0;
% phi2 = pi;
%     
% each_pre_t = 5*1000;
% each_stim_t = 6*1000; 
% each_break_t = 4*1000;
% total_t = 61*1000; 
% 
% stim_tt = dt:dt:each_stim_t;
% 
% I1_stim(a,:) = A1*cos(2*pi*f1*stim_tt+phi1);
% I2_stim(a,:) = A2*cos(2*pi*f2*stim_tt+phi1);
% 
% if ramp_up_t
%     ramp_up_size = round(ramp_up_t/dt);
%     ramp_vec = 0:1/ramp_up_size:1;
% 
%     I1_stim(a,1:length(ramp_vec)) = ramp_vec.*I1_stim(a,1:length(ramp_vec));
%     I2_stim(a,1:length(ramp_vec)) = ramp_vec.*I2_stim(a,1:length(ramp_vec));
% end
% 
% if ramp_down_t
%     ramp_down_size = round(ramp_down_t/dt);
%     ramp_vec = fliplr(0:1/ramp_down_size:1);
% 
%     I1_stim(a,end-length(ramp_vec)+1:end) = ramp_vec.*I1_stim(a,end-length(ramp_vec)+1:end);
%     I2_stim(a,end-length(ramp_vec)+1:end) = ramp_vec.*I2_stim(a,end-length(ramp_vec)+1:end);
%     
% end
% 
% % figure,plot(I1_stim(a,:)+ I2_stim(a,:))
% 
% % period of zeros
% 
% break_tt = dt:dt:each_break_t;
% 
% I1_break= zeros(1,length(break_tt));
% I2_break= zeros(1,length(break_tt));
% 
% % pre-stim period
% 
% pre_tt = dt:dt:each_pre_t;
% 
% I1_pre = zeros(1,length(pre_tt));
% I2_pre = zeros(1,length(pre_tt));
% 
% % complete 'cycle'
% 
% I1_cycle(a,:) = [I1_pre  I1_stim(a,:), I1_break];
% 
% I2_cycle(a,:) = [I2_pre  I2_stim(a,:), I2_break];
% 
% % figure,plot(I1_cycle.(char(sf(a))))
% 
% I1x = [I1x, I1_cycle(a,:)];
% I2x = [I2x, I2_cycle(a,:)];
% end
% 
% figure,plot(I1x + I2x)
% % saveas(gcf,'quantum_1ep', '.fig');
% 
% I1 = I1x;
% I2 = I2x;
% save('quantum_2ep','I1','I2');

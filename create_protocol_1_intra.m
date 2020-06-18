close all
clear all
clc

%% Stimulation signal script
df = 10;
sdf = mat2str(df);
carrier_f = [df, 7, 47, 97, 997, 2497, 4997];

for i = 1:length(carrier_f)
scf = mat2str(carrier_f(1,i));
if carrier_f(1,i) < 50
    amp(:,1) = [0.25:0.25:2]';
    amp(:,2) = [0.25:0.25:2]';
elseif carrier_f(1,i) > 50
    amp(:,1) = [1.5:0.05:2]';
    amp(:,2) = [1.5:0.05:2]';
end

random_amp_pairs_1= amp;
s = strcat('1_df_',sdf,'_1ep_',scf);
ss = strcat('1_df_',sdf,'_1ep_intra_',scf);

one_protocol_1 = random_amp_pairs_1;
save(s, 'one_protocol_1');

%% Create waveforms 
% 1 electrode pair - quantum

I1x = [];

dt = 0.004;
pulse_f = df;

f1 = carrier_f(1,i)/1000;
if carrier_f(1,i) == df
    f2 = df/1000;
else
    f2 = (carrier_f(1,i) + pulse_f)/1000;
end

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
save(ss,'I1','I2');
clearvars amp
end

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

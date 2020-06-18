close all
clear all
clc

%% Stimulation signal script

% Create stimulation protocol_1 with a random order 

freq_pairs = [10, 10; 50, 50; 100, 100; 500, 500; 1000, 1000; 2500, 2500; 5000, 5000];   
random_freq_pairs = freq_pairs(randperm(size(freq_pairs, 1)), :);

control = random_freq_pairs;

save('control', 'control');

%% Create waveforms

% Single-electrode pair - protocol 1

I1 = [];

for i = 1:length(freq_pairs)
        
        dt = 0.004;
        f1(i,1) = random_freq_pairs(i,1);
        f2(i,1) = random_freq_pairs(i,2);
        
        ramp_up_t = 500; % ms
        ramp_down_t = 500;
        
        A1 = 0.5;
        A2 = 0.5;
        
        phi1 = 0;
        phi2 = pi;
        
        % convert to ms
        
        f1(i,1) = f1(i,1)/1000;
        f2(i,1) = f2(i,1)/1000;
        
        each_stim_t = 6*1000; % 10 s on
        each_break_t = 5*1000;
        total_t = 90*1000; % 1 min
        
        stim_tt = dt:dt:each_stim_t;
        
        I1_stim(i,:) = A1*cos(2*pi*f1(i,1)*stim_tt+phi1) + A2*cos(2*pi*f2(i,1)*stim_tt+phi1);
        
        if ramp_up_t
            ramp_up_size = round(ramp_up_t/dt);
            ramp_vec = 0:1/ramp_up_size:1;
            
             I1_stim(i,1:length(ramp_vec)) = ramp_vec.*I1_stim(i,1:length(ramp_vec));
            
        end
        
        if ramp_down_t
            ramp_down_size = round(ramp_down_t/dt);
            ramp_vec = fliplr(0:1/ramp_down_size:1);
            
             I1_stim(i,end-length(ramp_vec)+1:end) = ramp_vec.*I1_stim(i,end-length(ramp_vec)+1:end);
            
        end
        
%         figure,plot( I1_stim(i,:))
        
        % period 10s of zeros
        break_tt = dt:dt:each_break_t;
        
        I1_break = zeros(1,length(break_tt));
        
        % complete 'cycle'
        I1_cycle(i,:) = [I1_stim(i,:) I1_break];
        
%         figure,plot(I1_cycle(i,:))
        
        I1 = [I1, I1_cycle(i,:)];

end
each_pre_t = 5*1000;
pre_tt = dt:dt:each_pre_t;
I1_pre = zeros(1,length(pre_tt));

I1 = [I1_pre, I1];
I2 = zeros(size(I1));
figure,plot(I1 + I2)
save('control_1ep' ,'I1','I2')
saveas(gcf,['control_1ep' '.fig']);

% %% Two electrode pair - protocol 1
% 
% I1 = [];
% I2 = [];
% 
% for i = 1:length(freq_pairs)
%         
%         dt = 0.004;
%         f1(i,1) = random_freq_pairs(i,1);
%         f2(i,1) = random_freq_pairs(i,2);
%         
%         ramp_up_t = 500; % ms
%         ramp_down_t = 500;
%         
%         A1 = 0.5;
%         A2 = 0.5;
%         
%         phi1 = 0;
%         phi2 = pi;
%         
%         % convert to ms
%         
%         f1(i,1) = f1(i,1)/1000;
%         f2(i,1) = f2(i,1)/1000;
%         
%         each_pre_t = 1*1000;
%         each_stim_t = 5*1000; % 10 s on
%         each_break_t = 4*1000;
%         total_t = 100*1000; % 1 min
%         
%         stim_tt = dt:dt:each_stim_t;
%         
%         I1_stim(i,:) = A1*cos(2*pi*f1(i,1)*stim_tt+phi1);
%         I2_stim(i,:) = A2*cos(2*pi*f2(i,1)*stim_tt+phi1);
%         
%         if ramp_up_t
%             ramp_up_size = round(ramp_up_t/dt);
%             ramp_vec = 0:1/ramp_up_size:1;
%             
%              I1_stim(1:length(ramp_vec)) = ramp_vec.*I1_stim(1:length(ramp_vec));
%              I2_stim(1:length(ramp_vec)) = ramp_vec.*I2_stim(1:length(ramp_vec));
%         end
%         
%         if ramp_down_t
%             ramp_down_size = round(ramp_down_t/dt);
%             ramp_vec = fliplr(0:1/ramp_down_size:1);
%             
%              I1_stim(end-length(ramp_vec)+1:end) = ramp_vec.*I1_stim(end-length(ramp_vec)+1:end);
%              I2_stim(end-length(ramp_vec)+1:end) = ramp_vec.*I2_stim(end-length(ramp_vec)+1:end);
%         end
%         
% %         figure,plot( I1_stim(i,:))
%         
%         % period 10s of zeros
%         break_tt = dt:dt:each_break_t;
%         
%         I1_break = zeros(1,length(break_tt));
%         I2_break = zeros(1,length(break_tt));
%         
%         % pre-stim period
%         
%         pre_tt = dt:dt:each_pre_t;
%         
%         I1_pre = zeros(1,length(pre_tt));
%         I2_pre = zeros(1,length(pre_tt));
%         
%         % complete 'cycle'
%         I1_cycle(i,:) = [I1_pre  I1_stim(i,:) I1_break];
%         I2_cycle(i,:) = [I2_pre  I2_stim(i,:) I2_break];
%         
% %         figure,plot(I1_cycle(i,:))
%         
%         I1 = [I1, I1_cycle(i,:)];
%         I2 = [I2, I2_cycle(i,:)];
% end
% 
% figure,plot(I1 + I2)
% save('control_2ep' ,'I1','I2')
% saveas(gcf,['control_2ep' '.fig']);
% 

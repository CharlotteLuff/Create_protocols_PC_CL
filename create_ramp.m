close all
clear all
clc

%% Stimulation signal script
% Use a subset of amplitudes
% Create stimulation protocol_1 with a random order 

freq_pairs = [10, 10; 50, 50; 100, 100; 500, 500; 1000, 1000; 2500, 2500; 5000, 5000];
random_freq_pairs = freq_pairs(randperm(size(freq_pairs, 1)), :);

ramp = random_freq_pairs; 

save('ramp', 'ramp');

%% Create waveforms

% Single-electrode pair - protocol 1

I1_n = [];

ramp_up_time = [500 1000 5000 10000 20000]; % ms

for i = 1:length(freq_pairs)
    for x = 1:length(ramp_up_time)
        
        dt = 0.004;
        f1(i,1) = random_freq_pairs(i,1);
        f2(i,1) = random_freq_pairs(i,2);
        
        ramp_up_t = ramp_up_time(1,x);
        ramp_down_t = 500;
        
        A1 = 0.5;
        A2 = 0.5;
        
        phI1_n = 0;
        phI2_n = pi;
        
        % convert to ms
        
        f1(i,1) = f1(i,1)/1000;
        f2(i,1) = f2(i,1)/1000;
        
        each_stim_t = ramp_up_t + (5*1000);
        each_break_t = 5*1000;
        total_t = 90*1000; % 1 min
        
        stim_tt = dt:dt:each_stim_t;
        
        I1_n_stim{1,x}(i,:) = A1*cos(2*pi*f1(i,1)*stim_tt+phI1_n) + A2*cos(2*pi*f2(i,1)*stim_tt+phI1_n);
        
        if ramp_up_t
            ramp_up_size = round(ramp_up_t/dt);
            ramp_vec = 0:1/ramp_up_size:1;
            
            I1_n_stim{1,x}(i,1:length(ramp_vec)) = ramp_vec.*I1_n_stim{1,x}(i,1:length(ramp_vec));
            
        end
        
        if ramp_down_t
            ramp_down_size = round(ramp_down_t/dt);
            ramp_vec = fliplr(0:1/ramp_down_size:1);
            
            I1_n_stim{1,x}(i,end-length(ramp_vec)+1:end) = ramp_vec.*I1_n_stim{1,x}(i,end-length(ramp_vec)+1:end);
            
        end
        
        %         figure,plot( I1_n_stim(i,:))
        
        % period 10s of zeros
        break_tt = dt:dt:each_break_t;
        
        I1_n_break = zeros(1,length(break_tt));
        
        % complete 'cycle'
        I1_n_cycle{1,x}(i,:) = [I1_n_stim{1,x}(i,:) I1_n_break];
        
        %         figure,plot(I1_n_cycle(i,:))
    
    end
end     

I1_n = horzcat(I1_n_cycle{1,:});
each_pre_t = 5*1000;
pre_tt = dt:dt:each_pre_t;
I1_n_pre = zeros(1,length(pre_tt));

for i = 1:size(I1_n,1)
    I1_n_new(i,:) = [I1_n_pre, I1_n(i,:)];
end
I1_n = I1_n_new;
I2_n = zeros(size(I1_n_new));

for i = 1:size(I1_n,1)
    figure
    plot(I1_n_new(i,:)+I2_n(i,:))
    I1 = I1_n_new(i,:);
    I2 = I2_n(i,:);
    s2 = mat2str(freq_pairs(i,1));
    s1 = 'ramp_';
    s3 = strcat(s1,s2);
    save(s3 ,'I1','I2')
    saveas(gcf,[s3 '.fig']);
end

    % %% Two electrode pair - protocol 1
    %
    % I1_n = [];
    % I2_n = [];
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
    %         phI1_n = 0;
    %         phI2_n = pi;
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
    %         I1_n_stim(i,:) = A1*cos(2*pi*f1(i,1)*stim_tt+phI1_n);
    %         I2_n_stim(i,:) = A2*cos(2*pi*f2(i,1)*stim_tt+phI1_n);
    %
    %         if ramp_up_t
    %             ramp_up_size = round(ramp_up_t/dt);
    %             ramp_vec = 0:1/ramp_up_size:1;
    %
    %              I1_n_stim(1:length(ramp_vec)) = ramp_vec.*I1_n_stim(1:length(ramp_vec));
    %              I2_n_stim(1:length(ramp_vec)) = ramp_vec.*I2_n_stim(1:length(ramp_vec));
    %         end
    %
    %         if ramp_down_t
    %             ramp_down_size = round(ramp_down_t/dt);
    %             ramp_vec = fliplr(0:1/ramp_down_size:1);
    %
    %              I1_n_stim(end-length(ramp_vec)+1:end) = ramp_vec.*I1_n_stim(end-length(ramp_vec)+1:end);
    %              I2_n_stim(end-length(ramp_vec)+1:end) = ramp_vec.*I2_n_stim(end-length(ramp_vec)+1:end);
    %         end
    %
    % %         figure,plot( I1_n_stim(i,:))
    %
    %         % period 10s of zeros
    %         break_tt = dt:dt:each_break_t;
    %
    %         I1_n_break = zeros(1,length(break_tt));
    %         I2_n_break = zeros(1,length(break_tt));
    %
    %         % pre-stim period
    %
    %         pre_tt = dt:dt:each_pre_t;
    %
    %         I1_n_pre = zeros(1,length(pre_tt));
    %         I2_n_pre = zeros(1,length(pre_tt));
    %
    %         % complete 'cycle'
    %         I1_n_cycle(i,:) = [I1_n_pre  I1_n_stim(i,:) I1_n_break];
    %         I2_n_cycle(i,:) = [I2_n_pre  I2_n_stim(i,:) I2_n_break];
    %
    % %         figure,plot(I1_n_cycle(i,:))
    %
    %         I1_n = [I1_n, I1_n_cycle(i,:)];
    %         I2_n = [I2_n, I2_n_cycle(i,:)];
    % end
    %
    % figure,plot(I1_n + I2_n)
    % save('control_2ep' ,'I1_n','I2_n')
    % saveas(gcf,['control_2ep' '.fig']);
    %

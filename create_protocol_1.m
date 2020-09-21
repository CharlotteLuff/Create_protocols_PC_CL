close all
clear all
clc

% User input

df = [5 10 40]; % define df
sf = {'df5' ,'df10','df40'}; % define df
s1 = {'1_df5_1ep' ,'1_df10_1ep','1_df40_1ep'}; % define df
s2 = {'1_df5_2ep' ,'1_df10_2ep','1_df40_2ep'}; % define df

%% Stimulation signal script

% Create stimulation protocol_1 with a random order 

for x = 1:length(df)
freq_pairs = [df(x), df(x); 7 , 7 + df(x); 47, 47 + df(x); 97, 97 + df(x); 497, 497 + df(x); 997, 997 + df(x); 2497, 2497 + df(x); 4997, 4997 + df(x)];   
random_freq_pairs_1.(char(sf(x))) = freq_pairs(randperm(size(freq_pairs, 1)), :);
% random_freq_pairs_1.(char(sf(x))) =[df(x), df(x); 7 , 7 + df(x); 47, 47 + df(x); 97, 97 + df(x); 497, 497 + df(x); 997, 997 + df(x); 2497, 2497 + df(x); 4997, 4997 + df(x)];
end

one_five_hz = random_freq_pairs_1.df5;
one_ten_hz = random_freq_pairs_1.df10;
one_forty_hz = random_freq_pairs_1.df40;

save('1_5hz', 'one_five_hz');
save('1_10hz', 'one_ten_hz');
save('1_40hz', 'one_forty_hz');

%% Create waveforms

% 1 electrode pair - protocol 1

for  x = 1:length(df)
I1x.(char(sf(x))) = [];
end

for i = 1:length(freq_pairs)
    for x = 1:length(df)
        
        dt = 0.004;
        f1.(char(sf(x)))(i,1) = random_freq_pairs_1.(char(sf(x)))(i,1);
        f2.(char(sf(x)))(i,1) = random_freq_pairs_1.(char(sf(x)))(i,2);
        
        ramp_up_t = 500; % ms
        ramp_down_t = 500;
        
        A1 = 0.5;
        A2 = 0.5;
        
        phi1 = 0;
        phi2 = pi;
        
        % convert to ms
        
        f1.(char(sf(x)))(i,1) = f1.(char(sf(x)))(i,1)/1000;
        f2.(char(sf(x)))(i,1) = f2.(char(sf(x)))(i,1)/1000;
        
        each_stim_t = 6*1000;
        each_break_t = 5*1000;
        total_t = 93*1000; 
        
        stim_tt = dt:dt:each_stim_t;
        
        I1_stim.(char(sf(x)))(i,:) = A1*cos(2*pi*f1.(char(sf(x)))(i,1)*stim_tt+phi1) + A2*cos(2*pi*f2.(char(sf(x)))(i,1)*stim_tt+phi1);

        if ramp_up_t
            ramp_up_size = round(ramp_up_t/dt);
            ramp_vec = 0:1/ramp_up_size:1;
            
             I1_stim.(char(sf(x)))(i,1:length(ramp_vec)) = ramp_vec.*I1_stim.(char(sf(x)))(i,1:length(ramp_vec));

        end
        
        if ramp_down_t
            ramp_down_size = round(ramp_down_t/dt);
            ramp_vec = fliplr(0:1/ramp_down_size:1);
            
             I1_stim.(char(sf(x)))(i,end-length(ramp_vec)+1:end) = ramp_vec.*I1_stim.(char(sf(x)))(i,end-length(ramp_vec)+1:end);
            
        end
        
%         figure,plot( I1_stim.(char(sf(x)))(i,:))
%         
        % period of zeros
        
        break_tt = dt:dt:each_break_t;
        
        I1_break = zeros(1,length(break_tt));
        
        % complete 'cycle'
        
        I1_cycle.(char(sf(x)))(i,:) = [I1_stim.(char(sf(x)))(i,:) I1_break];
        
%         figure,plot(I1_cycle.(char(sf(x)))(i,:))
                % pre-stim period
                
        I1x.(char(sf(x))) = [I1x.(char(sf(x))), I1_cycle.(char(sf(x)))(i,:)];

    end
end

each_pre_t = 5*1000;
pre_tt = dt:dt:each_pre_t;
I1_pre = zeros(1,length(pre_tt));
        
for  x = 1:length(df)
    I1x.(char(sf(x))) = [I1_pre, I1x.(char(sf(x)))];
    I2x.(char(sf(x))) = zeros(size(I1x.(char(sf(x)))));
    figure,plot(I1x.(char(sf(x))) + I2x.(char(sf(x))))
    saveas(gcf,[char(s1(x)) '.fig']);
end

for  x = 1:length(df)
I1 = I1x.(char(sf(x)));
I2 = I2x.(char(sf(x)));
save(char(s1(x)),'I1','I2');
end
% 
% %% 2 electrode pair - protocol 1
% % fix unknown issue
% for  x = 1:length(df)
% I1x.(char(sf(x))) = [];
% I2x.(char(sf(x))) = [];
% end
% 
% for i = 1:length(freq_pairs)
%     for x = 1:length(df)
%         
%         dt = 0.004;
%         f1.(char(sf(x)))(i,1) = random_freq_pairs_1.(char(sf(x)))(i,1);
%         f2.(char(sf(x)))(i,1) = random_freq_pairs_1.(char(sf(x)))(i,2);
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
%         f1.(char(sf(x)))(i,1) = f1.(char(sf(x)))(i,1)/1000;
%         f2.(char(sf(x)))(i,1) = f2.(char(sf(x)))(i,1)/1000;
%         
%         each_pre_t = 1*1000;
%         each_stim_t = 6*1000; % 10 s on
%         each_break_t = 4*1000;
%         total_t = 101*1000; % 1 min
%         
%         stim_tt = dt:dt:each_stim_t;
%         
%         I1_stim.(char(sf(x)))(i,:) = A1*cos(2*pi*f1.(char(sf(x)))(i,1)*stim_tt+phi1);
%         I2_stim.(char(sf(x)))(i,:) = A2*cos(2*pi*f2.(char(sf(x)))(i,1)*stim_tt+phi1);
%         
%         if ramp_up_t
%             ramp_up_size = round(ramp_up_t/dt);
%             ramp_vec = 0:1/ramp_up_size:1;
%             
%              I1_stim.(char(sf(x)))(i,1:length(ramp_vec)) = ramp_vec.*I1_stim.(char(sf(x)))(i,1:length(ramp_vec));
%              I2_stim.(char(sf(x)))(i,1:length(ramp_vec)) = ramp_vec.*I2_stim.(char(sf(x)))(i,1:length(ramp_vec));
%         end
%         
%         if ramp_down_t
%             ramp_down_size = round(ramp_down_t/dt);
%             ramp_vec = fliplr(0:1/ramp_down_size:1);
%             
%              I1_stim.(char(sf(x)))(i,end-length(ramp_vec)+1:end) = ramp_vec.*I1_stim.(char(sf(x)))(i,end-length(ramp_vec)+1:end);
%              I2_stim.(char(sf(x)))(i,end-length(ramp_vec)+1:end) = ramp_vec.*I2_stim.(char(sf(x)))(i,end-length(ramp_vec)+1:end);
%         end
%         
% %         figure,plot(I1_stim.(char(sf(x)))(i,:) + I2_stim.(char(sf(x)))(i,:))
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
%         I1_cycle.(char(sf(x)))(i,:) = [I1_pre  I1_stim.(char(sf(x)))(i,:) I1_break];
%         I2_cycle.(char(sf(x)))(i,:) = [I2_pre  I2_stim.(char(sf(x)))(i,:) I2_break];
%         
% %         figure,plot(I1_cycle.(char(sf(x)))(i,:) + I2_cycle.(char(sf(x)))(i,:))
%         
%         I1x.(char(sf(x))) = [I1x.(char(sf(x))), I1_cycle.(char(sf(x)))(i,:)];
%         I2x.(char(sf(x))) = [I2x.(char(sf(x))), I2_cycle.(char(sf(x)))(i,:)];
%     end
% end
% 
% for  x = 1:length(df)
%     figure,plot(I1x.(char(sf(x))) + I2x.(char(sf(x))))
%     saveas(gcf,[char(s2(x)) '.fig']);
% end
% 
% for  x = 1:length(df)
% I1 = I1x.(char(sf(x)));
% I2 = I2x.(char(sf(x)));
% save(char(s2(x)),'I1','I2');
% end
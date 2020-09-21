close all;
clear;
clc;

%%%%%%%%%%%%%% TI parameters %%%%%%%%%%%%%%
f_TI_1 = 2000; 
f_TI_2 = 4000;                        %[editable-stim] 4000 for 0.5s; 2500 for 2ms; 2100 for 10ms
phi1 = 0;
phi2 = pi;

%%%%%%%%%%%%%% temporal parameters %%%%%%%%%%%%%%
dt= 0.004; %ms  0.004;  
pre_time = 1000; %ms                   [editable-time]
post_time = 1000; %ms                  [editable-time]
ramp_up = 500; %ms                     [editable-time]
ramp_down = 500; %ms                   [editable-time]

p_d =1000/abs(f_TI_1-f_TI_2);          %pulse_duration (ms) [uneditable-stim]
pulse_number_1=1;                      %[editable-stim]
pulse_number_4=4;
IPI=0; %inter-pulse interval (ms)      [editable-stim]

resting=4000; %ms resting time among sections [editable-time]

%%%%%%%%%%%%%% amplitude parameters %%%%%%%%%%%%%%
basis_amp= 0;                           %[editable-stim]
step_amp=0.25;                         %[editable-stim]
amp_scope = 8;                         %[editable-stim] e.g. from 0 to 2 mA, or from 2 to 4 mA
amp_num = amp_scope/step_amp;          %[uneditable-stim] 

%%%%%%%%%%%%%%  (1) pre stimulation  %%%%%%%%%%%%%%
l_p_s = pre_time/dt; % length points
pre_amp(1:1:l_p_s) = (1:1:l_p_s)*0;  %pre_time ampplitude

%%%%%%%%%%%%%%  (2) resting time     %%%%%%%%%%%%%%
l_p_s = resting/dt; % length points
rest_1_amp(1:1:l_p_s) = 0;

%%%%%%%%%%%%%%  (3) TI stimulation   %%%%%%%%%%%%%%
for j=1:1:amp_num % amplitude number    
    amp_scale=basis_amp+j*step_amp;
%%%  (3-1) TI 1 pulse   %%%
    %ramp-up
    l_p_s = ramp_up/dt; % length points
    ramp_up_amp_TI_1(1:1:l_p_s) = amp_scale*0.5*(sin(2*pi*f_TI_1*(1:1:l_p_s)*dt/1000+phi1) + sin(2*pi*f_TI_2*(1:1:l_p_s)*dt/1000+phi2));
    ramp_up_amp_TI_1 = ramp_up_amp_TI_1.*((1:1:l_p_s)/l_p_s); %*ramp_up_factor_TI
    
    %stim
    for i=1:1:pulse_number_1
        l_p_s = (p_d+IPI)/dt; % length points
        l_p_s_stim = p_d/dt; % length points
        stim_amp_TI_1( 1+(i-1)*l_p_s :1: l_p_s_stim+l_p_s*(i-1)) = amp_scale*0.5*(sin(2*pi*f_TI_1*(1:1:l_p_s_stim)*dt/1000+phi1) + sin(2*pi*f_TI_2*(1:1:l_p_s_stim)*dt/1000+phi2)); 
    end
    
    %ramp-down
    l_p_s = ramp_down/dt; % length points
    ramp_down_amp_TI_1(1:1:l_p_s) = amp_scale*0.5*(sin(2*pi*f_TI_1*(1:1:l_p_s)*dt/1000+phi1) + sin(2*pi*f_TI_2*(1:1:l_p_s)*dt/1000+phi2));
    ramp_down_amp_TI_1 = ramp_down_amp_TI_1.*((l_p_s:-1:1)/l_p_s); %*ramp_down_factor_TI
    
%%%  (3-2) TI 4 pulse   %%%%    
    %ramp-up
    l_p_s = ramp_up/dt; % length points
    ramp_up_amp_TI_4(1:1:l_p_s) = amp_scale*0.5*(sin(2*pi*f_TI_1*(1:1:l_p_s)*dt/1000+phi1) + sin(2*pi*f_TI_2*(1:1:l_p_s)*dt/1000+phi2));
    ramp_up_amp_TI_4 = ramp_up_amp_TI_4.*((1:1:l_p_s)/l_p_s); %*ramp_up_factor_TI
    
    %stim
    for i=1:1:pulse_number_4
        l_p_s = (p_d+IPI)/dt; % length points
        l_p_s_stim = p_d/dt; % length points
        stim_amp_TI_4( 1+(i-1)*l_p_s :1: l_p_s_stim+l_p_s*(i-1)) = amp_scale*0.5*(sin(2*pi*f_TI_1*(1:1:l_p_s_stim)*dt/1000+phi1) + sin(2*pi*f_TI_2*(1:1:l_p_s_stim)*dt/1000+phi2)); 
            if i ~= pulse_number_4
                stim_amp_TI_4( l_p_s_stim+(i-1)*l_p_s+1 :1: l_p_s+l_p_s*(i-1)) = 0; %used to set IPI
            end
    end
    
    %ramp-down
    l_p_s = ramp_down/dt; % length points
    ramp_down_amp_TI_4(1:1:l_p_s) = amp_scale*0.5*(sin(2*pi*f_TI_1*(1:1:l_p_s)*dt/1000+phi1) + sin(2*pi*f_TI_2*(1:1:l_p_s)*dt/1000+phi2));
    ramp_down_amp_TI_4 = ramp_down_amp_TI_4.*((l_p_s:-1:1)/l_p_s); %*ramp_down_factor_TI
    
    %all stim
    if j < amp_num % add resting time between two stimualtion sections
    stim_amp_TI_help = [ramp_up_amp_TI_1 stim_amp_TI_1 ramp_down_amp_TI_1 rest_1_amp ramp_up_amp_TI_4 stim_amp_TI_4 ramp_down_amp_TI_4 rest_1_amp];
    else 
        stim_amp_TI_help = [ramp_up_amp_TI_1 stim_amp_TI_1 ramp_down_amp_TI_1 rest_1_amp ramp_up_amp_TI_4 stim_amp_TI_4 ramp_down_amp_TI_4 ];
    end
    
    if 1 == j % combine all stimulation sections
        stim_amp_TI_array = [stim_amp_TI_help];
    else
        stim_amp_TI_array = [stim_amp_TI_array stim_amp_TI_help];
    end
end

%%%%%%%%%%%%%%  (4) post stimulation  %%%%%%%%%%%%%%
l_p_s = post_time/dt;
post_amp(1:1:l_p_s) = (1:1:l_p_s)*0;  %pre_time ampplitude

%%%%%%%%%%%%%%  Plot and store  %%%%%%%%%%%%%%
I1 = [pre_amp stim_amp_TI_array post_amp];
I2 = I1.*0;
all_pulse_stim_time=(1:1:length(I1))*dt/1000;

hold on;
plot(all_pulse_stim_time,I1,'k');
hold on;
plot(all_pulse_stim_time,I1,'.k');
xlabel('Time(s)');
ylabel('Amplitude');
save('all_pulse_stim_amp','I1','I2');
% axis([35 38.5 -1.1 1.1]);
%% Stimulation signal script
dt = 0.01;
pulse_f = [20, 40, 60, 80, 100];
ramp_up_t = 500; % ms
ramp_down_t = 500;

A1 = 0.5;
A2 = 0.5;

phi1 = 0;
phi2 = pi;
    
% convert to ms

for i = 1:length(pulse_f)
spulse_f = mat2str(pulse_f(1,i));
pulse_f(1,i) = pulse_f(1,i)/1000;

    f1 = pulse_f(1,i);
    f2 = pulse_f(1,i);
    
    phi1 = -pi/2;
    phi2 = -pi/2;
    
each_stim_t = 60*1000; % 10 s on
total_t = 60*1000; % 1 min

stim_tt = dt:dt:each_stim_t;

I1_stim = A1*cos(2*pi*f1*stim_tt+phi1);
I2_stim = A2*cos(2*pi*f2*stim_tt+phi2);

if ramp_up_t
    ramp_up_size = round(ramp_up_t/dt);
    ramp_vec = 0:1/ramp_up_size:1;

    I1_stim(1:length(ramp_vec)) = ramp_vec.*I1_stim(1:length(ramp_vec));
    I2_stim(1:length(ramp_vec)) = ramp_vec.*I2_stim(1:length(ramp_vec));

end

if ramp_down_t
    ramp_down_size = round(ramp_down_t/dt);
    ramp_vec = fliplr(0:1/ramp_down_size:1);

    I1_stim(end-length(ramp_vec)+1:end) = ramp_vec.*I1_stim(end-length(ramp_vec)+1:end);
    I2_stim(end-length(ramp_vec)+1:end) = ramp_vec.*I2_stim(end-length(ramp_vec)+1:end);
end
% 
% figure,plot(I1_stim+I2_stim)

% period 10s of zeros

% complete 'cycle'
I1_cycle = [I1_stim];
I2_cycle = [I2_stim];

% figure,plot(I1_cycle+I2_cycle)

%%
% full stim waveform
no_cycles = total_t/(each_stim_t);

I1 = repmat(I1_cycle,1,no_cycles);
I2 = repmat(I2_cycle,1,no_cycles);

% figure,plot(I1+I2)

%% save
s = strcat('microglia_waveform_1min','_1ep_',spulse_f,'hz.mat');
save(s,'I1','I2')

end


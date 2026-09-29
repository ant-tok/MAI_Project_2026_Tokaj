clear; close all; clc;
DAS_points = 500;
distance = 1.65;
c = 300000000;
data = readmatrix('ranging_data.csv');
x = data(:,1);
y = data(:,2:4);
time_diff = x(2) - x(1);
%d = x*c;
d = x.*c;
in_sig = data(:,2);
a1_sig = data(:,3);
a2_sig = data(:,4);


figure
plot(x, a1_sig, 'LineWidth', 1.2); hold on
plot(x, a2_sig, 'LineWidth', 1.2)
title('Raw antenna signals')
xlabel('Time (s)');
ylabel('Amplitude') 
legend('antenna 1', 'antenna 2', 'Location', 'northwest')
grid on

% Assumes peaks are caused from reflection. Not robust to noise (maybe
% needs smoothing)
[~, a1_idx] = max(a1_sig);
[~, a2_idx] = max(a2_sig);


pos_shift = 1;
if a1_idx < a2_idx
    disp('antenna 1 arrived first');
elseif a1_idx > a2_idx
    disp('antenna 2 arrived first');
    pos_shift = 0;
else
    disp('no measurable lead');
end



N = length(a1_sig);

X1 = fft(a1_sig, N); 
X2 = fft(a2_sig, N);   

%Get cross corelations
a1_auto_correlation = sum(a1_sig.^2);  
a2_auto_correlation = sum(a2_sig.^2);

C_Correlation = real(ifft(conj(X1) .* X2));


C_centered = fftshift(C_Correlation);
lags = (-floor(N/2):(ceil(N/2)-1))'; 



figure
plot(lags*time_diff*c/2, C_centered, 'LineWidth', 1.2)

title('Cross Correlation vs time delay')
xlabel('Distance difference')
ylabel('Cross-Correlation');
grid on

max_shift_samples = ceil((2*distance/c) / time_diff);
max_shift_samples = min(max_shift_samples, floor(N/2));
center = floor(N/2) + 1; 

if pos_shift
    win = center:(center + max_shift_samples - 1);
else
    win = (center - max_shift_samples):center;
end

e_array1_freq = a1_auto_correlation + a2_auto_correlation + 2*C_centered(win);
sig_dis_freq  = lags(win) * time_diff * c / 2;

[peak_energy_f, peak_idx_f] = max(e_array1_freq);
best_distance_f = sig_dis_freq(peak_idx_f);
best_k_f = lags(win(peak_idx_f));



% figure
% plot(sig_dis_freq, e_array1_freq, 'LineWidth', 1.2)
% hold on
% 
% title('Summed energy vs distance ')
% xlabel('Relative distance')
% ylabel('Energy')
% grid on


% Get distance from antenna (assumes reflector directly in between)
distance_from_antenna = (distance - best_distance_f)/2

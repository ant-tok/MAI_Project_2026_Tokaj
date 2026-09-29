clear; close all; clc;
DAS_points = 50;
distance = 1.65;
c = 300000000;
data = readmatrix('ranging_data.csv');
x = data(:,1);
y = data(:,2:4);
time_diff = x(2) - x(1);
d = x.*c;
in_sig = data(:,2);
a1_sig = data(:,3);
a2_sig = data(:,4);


figure
plot(x, a1_sig, 'LineWidth', 1.2); hold on
plot(x, a2_sig, 'LineWidth', 1.2)
title('Raw antenna signals (before alignment)')
xlabel('Time (s)');
ylabel('Amplitude')
legend('antenna 1', 'antenna 2', 'Location', 'northwest')
grid on

[~, a1_idx] = max(a1_sig);
[~, a2_idx] = max(a2_sig);

pos_shift = 1;
if a1_idx < a2_idx
    lead_msg_f = 'antenna 1 arrived first';
elseif a1_idx > a2_idx
    lead_msg_f = 'antenna 2 arrived first';
    pos_shift = 0;
else
    lead_msg_f = 'no measurable lead';
end


N = length(a1_sig);
X1 = fft(a1_sig, N);    
m = (0:N-1)';      

E2 = sum(a2_sig.^2);    

max_shift_samples = ceil((2*distance/c) / time_diff);
max_shift_samples = min(max_shift_samples, floor(N/2));

if pos_shift
    k_range = 0:(max_shift_samples - 1);
else
    k_range = -max_shift_samples:0;
end

n_k = length(k_range);
e_array1_freq = zeros(n_k, 1);
sig_dis_freq  = zeros(n_k, 1);

for idx = 1:n_k
    k = k_range(idx);


    delta_spectrum = exp(-1i * 2*pi * k * m / N);
    shifted_a1 = real(ifft(X1 .* delta_spectrum));

    e_array1_freq(idx) = sum((shifted_a1 + a2_sig).^2);
    sig_dis_freq(idx)  = k * time_diff * c / 2;
end

[peak_energy_f, peak_idx_f] = max(e_array1_freq);
best_distance_f = sig_dis_freq(peak_idx_f);
best_k_f = k_range(peak_idx_f);


figure
plot(sig_dis_freq, e_array1_freq, 'LineWidth', 1.2)
hold on

xlabel('Relative distance, a1 vs a2 (signed, m)')
ylabel('Energy')

grid on


% Get distance from antenna
distance_from_antenna = (distance - best_distance_f)/2
clear; close all; clc;

DAS_points = 500;

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
plot(d,y, 'LineWidth', 1.2) 
xlabel('Distance')
ylabel('Signal')
legend('In signal','Ant 1','Ant 2','Location','southeast')
grid on

% compute FFT and frequency axis
Fs = 1 / time_diff;            % sampling frequency
N = length(in_sig);
nfft = 2^nextpow2(N);
f = (0:(nfft/2))' * (Fs / nfft);

% FFTs
IN = fft(in_sig, nfft);
A1 = fft(a1_sig, nfft);
A2 = fft(a2_sig, nfft);

% single-sided magnitude (normalized) and phase (unwrap)
mag_in = abs(IN(1:nfft/2+1)) / N;
mag_a1 = abs(A1(1:nfft/2+1)) / N;
mag_a2 = abs(A2(1:nfft/2+1)) / N;

phase_in = unwrap(angle(IN(1:nfft/2+1)));
phase_a1 = unwrap(angle(A1(1:nfft/2+1)));
phase_a2 = unwrap(angle(A2(1:nfft/2+1)));


% plot magnitude
figure
plot(f(1:200), mag_in(1:200), 'LineWidth',1.2); hold on
plot(f(1:200), mag_a1(1:200), 'LineWidth',1.2)
plot(f(1:200), mag_a2(1:200), 'LineWidth',1.2)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('in\_sig','a1\_sig','a2\_sig','Location','best')
grid on

% plot phase
figure
plot(f(1:500), phase_in(1:500), 'LineWidth',1.2); hold on
plot(f(1:500), phase_a1(1:500), 'LineWidth',1.2)
plot(f(1:500), phase_a2(1:500), 'LineWidth',1.2)
xlabel('Frequency (Hz)')
ylabel('Phase  (rads)')
legend('in\_sig','a1\_sig','a2\_sig','Location','best')
grid on

samples = length(a1_sig);

a = in_sig + a1_sig;


squared = a .^ 2;

e = sum(squared);


b = samples / DAS_points;
shift_time = time_diff * b;
e_array1 = zeros(DAS_points,1);
e_array2 = zeros(DAS_points,1);
for i = 1:DAS_points
    
    k = round(i*b);
    %shift_sig1 = zeros(size(a1_sig));
    %shift_sig2 = zeros(size(a2_sig));

    shift_sig1 = circshift(a1_sig, -k);
    shift_sig2 = circshift(a2_sig, -k);

    

    summed1 = shift_sig1 + in_sig;
    squared1 = summed1 .^ 2;

    summed2 = shift_sig2 + in_sig;
    squared2 = summed2 .^ 2;

    e1 = sum(squared1);
    e2 = sum(squared2);
    
    
    e_array1(i) = e1;
    e_array2(i) = e2;
end




step = (shift_time * c) / 2;

to_distance = round(distance/step) + 1;


summed_energy = zeros(to_distance,1);



n1 = min(DAS_points, to_distance);
n2 = max(DAS_points, to_distance);

sig_dis = (0:(n2-1))' * (shift_time * c) / 2;

summed_energy(1:n1) = e_array1(1:n1);

% remaining points to fill from reversed e_array2
rem = to_distance - n1;
if rem > 0
    rev_e2 = flipud(e_array2);
    summed_energy(n2-n1+1:n2) = rev_e2(1:n1);
    
end





figure
plot(sig_dis,summed_energy, 'LineWidth', 1.2)  

title('Summed Energy')

xlabel('Distance') 
ylabel('Energy');
grid on

ave_energy = mean(summed_energy);
abs_energy = abs(summed_energy - ave_energy);

figure
plot(sig_dis,abs_energy, 'LineWidth', 1.2)  

title('Distance from Mean Energy')
xlabel('Distance') 
ylabel('Energy');
grid on


window_len = 20;
sm_abs_energy = movmean(abs_energy, window_len);


figure
plot(sig_dis, sm_abs_energy, 'LineWidth',1.2)
title('Smoothed Distance from Mean Energy')
xlabel('Distance')
ylabel('Energy')
grid on


heatmat = repmat(summed_energy(:),1,DAS_points);

figure
imagesc(sig_dis, 1:1, heatmat')
set(gca,'YDir','normal')
xlabel('Distance') 

title('Summed Energy Heatmap')
colorbar
colormap(jet)


heatmat = repmat(abs_energy(:),1,DAS_points);

figure
imagesc(sig_dis, 1:1, heatmat')
set(gca,'YDir','normal')
xlabel('Distance') 

title('Distance from Mean Energy Heatmap')
colorbar
colormap(jet)


heatmat = repmat(sm_abs_energy(:),1,DAS_points);

figure
imagesc(sig_dis, 1:1, heatmat')
set(gca,'YDir','normal')
xlabel('Distance') 

title('Distance from Mean Energy Heatmap')
colorbar
colormap(jet)
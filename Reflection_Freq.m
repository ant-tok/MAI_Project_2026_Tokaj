clear; close all; clc;
%model arm, look at TOF

fs_permittivity = 8.854187817e-12; % permittivity of free space 
fs_permeability = pi*4e-7; % permeability of free space 


rel_permittivity_array = [1, 4, 9, 2];   
layer_thickness = [1.15, 1.5, 2.0]; 

N = numel(rel_permittivity_array);

dist_sum = zeros(1, N-1);
for n = 1:N-1
    dist_sum(n) = sum(layer_thickness(1:n));
end


zmax = 10.0;
tim_max = 80;

frequency_spacing = 1e7;
max_freq = 3e9; 
no_tones = max_freq/frequency_spacing;
freq_list = (1:no_tones)'*frequency_spacing;

% Values for pulse
f0 = 1e9; %Carrier frequency
tp = 1e-9; % Envelope width


%c0 = 1/sqrt(fs_permittivity*fs_permeability); 
c0 = 3e8;


v_array = c0 ./ sqrt(rel_permittivity_array);
imped_array = sqrt(fs_permeability ./ (fs_permittivity*rel_permittivity_array));
 

Gamma_array = zeros(1, N-1);
Tau_array = zeros(1, N-1);
for m = 1:N-1
    eta_l = imped_array(m);
    eta_r = imped_array(m+1);
    Gamma_array(m) = (eta_r - eta_l)/(eta_r + eta_l);
    Tau_array(m) = 2*eta_r/(eta_r + eta_l);
end

tau_total = ones(1, N);
for n = 2:N
    tau_total(n) = tau_total(n-1)*Tau_array(n-1);
end

no_tones = max_freq/frequency_spacing;

Nf = numel(freq_list);
w = 2*pi*freq_list;
 
k_array = w * (1./v_array);
k1 = k_array(:,1);
 

 
phase_at_bourdary = zeros(Nf, N-1);   % phi_to_boundary(:,m) = one-way phase to reach boundary m
running_phase = zeros(Nf,1);
for m = 1:N-1
    running_phase = running_phase + k_array(:,m).*layer_thickness(m);
    phase_at_bourdary(:,m) = running_phase;
end

dist_step = (v_array(1)/f0)/40;
z = (0:dist_step:zmax);
Nz = numel(z);
 

E_fwd = zeros(Nf, Nz);

 
for n = 1:N
    kn = k_array(:,n);
    if n == 1
        mask = z < dist_sum(1);
        E_fwd(:,mask) = tau_total(n) * exp(-1j*kn*(z(mask)));
    else
        if n == N
            mask = z >= dist_sum(N-1);
            z0 = dist_sum(N-1);
        else
            mask = (z >= dist_sum(n-1)) & (z < dist_sum(n));
            z0 = dist_sum(n-1);
        end
        E_fwd(:,mask) = tau_total(n) * exp(-1j*phase_at_bourdary(:,n-1)) .* exp(-1j*kn*(z(mask)-z0));

    end

    
end


E_echo = zeros(Nf, Nz);
for m = 1:N-1
    echo_amp = tau_total(m) * Gamma_array(m);   % refection after transmission
 
    for n = 1:m   % every region the echo must cross on its way back to the receiver
        if n == 1
            mask = z < dist_sum(1);
        else
            mask = (z >= dist_sum(n-1)) & (z < dist_sum(n));
        end

        a = echo_amp * exp(-1j*(2*phase_at_bourdary(:,m) - phase_at_bourdary(:,n)));
        a = a .* exp(-1j*k_array(:,n).*(dist_sum(n) - z(mask)));
        E_echo(:,mask) = E_echo(:,mask) +  a;
        % E_echo(:,mask) = E_echo(:,mask) + echo_amp * exp(-1j*phase_total);
    end
end
 
E = E_fwd + E_echo;


Nz = numel(z);

Y = E(:, 1);


Frequency_Weights = (tp*sqrt(pi)/2) * exp(-1j*2*pi*freq_list*3*tp);

Frequency_Weights = Frequency_Weights .* exp(-(pi*(freq_list - f0)*tp).^2);   % Gaussian bell around f0

response = Frequency_Weights.*Y;

Ifft_Length = 2^12;

time_axis = (0:Ifft_Length-1)'/(Ifft_Length*frequency_spacing);  


synth = @(S) real(ifft([0; Ifft_Length*frequency_spacing*S; zeros(Ifft_Length-2*no_tones-1,1); conj(flipud(Ifft_Length*frequency_spacing*S))]));
full_synth = @(S) real(ifft([zeros(1,Nz); Ifft_Length*frequency_spacing*S; zeros(Ifft_Length-2*no_tones-1, Nz);conj(flipud(Ifft_Length*frequency_spacing*S))]));

recieved = synth(Y);
transmitted = synth(Frequency_Weights); 
response_signal = synth(response);

response_full = Frequency_Weights .* E; 
Etz = full_synth((response_full));


% figure;
% 
% plot(freq_list/1e9, abs(Frequency_Weights), 'LineWidth', 1.2); grid on
% xlabel('Frequency (GHz)'); ylabel('|weight|');
% title('Pulse Spectrum');
% 
% 
% figure;
subplot(3,1,1);

plot(time_axis*1e9, recieved); grid on
xlim([0 tim_max]); xlabel('t (ns)'); ylabel('E (V/m)');
title('Recieved Signal in Time');

subplot(3,1,2);

plot(time_axis*1e9, transmitted); grid on
xlim([0 tim_max]); xlabel('t (ns)'); ylabel('E (V/m)');
title('Pulse');

subplot(3,1,3);

plot(time_axis*1e9, response_signal); grid on
xlim([0 tim_max]); xlabel('t (ns)'); ylabel('E at Rx (V/m)');
title('Recieved Conjugated with pulse');


figure;
imagesc(z, time_axis*1e9, Etz);
axis xy;
colormap(jet);
colorbar;
ylim([0 tim_max]); xlabel('distance (m)'); ylabel('time (ns)');
title('E at distance/time');

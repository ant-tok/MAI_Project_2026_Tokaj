clear; close all; clc;
%model arm, look at TOF
% It'is database,

function E = recursiveField(n, mask_area, k_array, dist_sum, imped_array)

    if n == 1
        E_entry = 1; 
        dist = 0;
    else
        E_prev_at_boundary = recursiveField(n-1, dist_sum(n-1), k_array, dist_sum, imped_array);

        tau = 2*imped_array(:,n) ./ (imped_array(:,n) + imped_array(:,n-1));

        E_entry = tau .* E_prev_at_boundary;
        dist = dist_sum(n-1);
    end

    E = E_entry .* exp(-1j*k_array(:,n).*(mask_area-dist));
end

function E = reflect_Multiplier(n, mask_area, m, Em, k_array, dist_sum, imped_array)
    % fprintf('%f\n', dist_sum(n));
    if n == m
        % fprintf("G\n")
        E_boundary = Em; 
    else
        tau = 2*imped_array(:,n) ./ (imped_array(:,n+1) + imped_array(:,n));
        E_boundary = tau .* reflect_Multiplier(n+1, dist_sum(n), m, Em, k_array, dist_sum, imped_array);
    end
    E = E_boundary .* exp(-1j*k_array(:,n).*(dist_sum(n) - mask_area));
end

function eps_r = debye_func(w, e_inf, e_s, Sig, fs_perm)

    eps_r = e_inf + (e_s - e_inf)./(1 + 1j*w*17.5*1e-12) - 1j*Sig./(w*fs_perm);
end
 


fs_permittivity = 8.854187817e-12; % permittivity of free space 
fs_permeability = pi*4e-7; % permeability of free space 

%Layer properties
% rel_permittivity_array = [2+1j, 9+1j, 1+1j, 4+1j];  
rel_permittivity_array = [2, 9, 1, 4];  



layer_thickness = [0.002, 0.05, 0.002]; %Treat last layer as infinite   

N = numel(rel_permittivity_array);

dist_sum = zeros(1, N-1);
for n = 1:N-1
    dist_sum(n) = sum(layer_thickness(1:n));
end


zmax = 0.1;
tim_max = 10;

frequency_spacing = 1e7;
max_freq = 3e9; 
no_Tones = max_freq/frequency_spacing;
freq_list = (1:no_Tones)'*frequency_spacing;


Nf = numel(freq_list);
w = 2*pi*freq_list;

two_d_permittivity_array = zeros(Nf,N);


skin_perm = debye_func(w, 18.51, 42.28, 0.679, fs_permittivity);
muscle_perm= debye_func(w, 24.60, 58.22, 0.664, fs_permittivity);

% two_d_permittivity_array(:,1) = 2;
two_d_permittivity_array(:,1) = skin_perm;
two_d_permittivity_array(:,2) = muscle_perm;
% two_d_permittivity_array(:,2) = 9;
% two_d_permittivity_array(:,3) = 1;
two_d_permittivity_array(:,3) = skin_perm;
two_d_permittivity_array(:,4) = 1;


%c0 = 1/sqrt(fs_permittivity*fs_permeability); 
c0 = 3e8;


imped_array = sqrt(fs_permeability ./ (fs_permittivity*rel_permittivity_array));


% Values for pulse
f0 = 1e9; 
tp = 25e-11;

dist_step = (c0/f0)/400;
z = (0:dist_step:zmax);
Nz = numel(z);


no_Tones = max_freq/frequency_spacing;

% v_array = c0 ./ sqrt(rel_permittivity_array);
% k_array = w * (1./v_array);


k_array = w .* sqrt(fs_permeability*fs_permittivity*two_d_permittivity_array);
imped_array = sqrt(fs_permeability ./ (fs_permittivity*two_d_permittivity_array));


% Find forward E, from transission coefficient
% Use of 8.88 from book and surrounding equations
E_fwd = zeros(Nf, Nz); 

for n = 1:N

    kn = k_array(:,n);

    if n == 1
        mask = z < dist_sum(1);
        E_fwd(:,mask) = recursiveField(n, z(mask), k_array, dist_sum, imped_array);

    else
        if n == N
          mask = z >= dist_sum(N-1);

        else
          mask = (z >= dist_sum(n-1)) & (z < dist_sum(n));

        end
        E_fwd(:,mask) = recursiveField(n, z(mask), k_array, dist_sum, imped_array);
    end
    
end


%Reflective E, or Backwards E
E_Back = zeros(Nf, Nz);
for m = 1:N-1


    leftside = imped_array(:,m);
    rightside = imped_array(:,m+1);
    Gamma_m = (rightside - leftside)./(rightside + leftside);
 
    E_at_m = Gamma_m .* recursiveField(m, dist_sum(m), k_array, dist_sum, imped_array);
 
    for n = 1:m 
        %fprintf("test1");
        if n == 1
            
            mask = z < dist_sum(1);
        else
            mask = (z >= dist_sum(n-1)) & (z < dist_sum(n));
        end
        edge_val = reflect_Multiplier(n, z(mask), m, E_at_m, k_array, dist_sum, imped_array);
        kn = k_array(:,n);
        E_Back(:,mask) = E_Back(:,mask) + edge_val;
    end
end

 

E = E_fwd + E_Back;


Nz = numel(z);

Y = E(:, 1);


Frequency_Weights = (tp*sqrt(pi)/2) * exp(-1j*2*pi*freq_list*6*tp);

Frequency_Weights = Frequency_Weights .* exp(-(pi*(freq_list - f0)*tp).^2);   % Gaussian bell around f0

response = Frequency_Weights.*Y;

Ifft_Length = 2^12;

time_axis = (0:Ifft_Length-1)'/(Ifft_Length*frequency_spacing);  


synth = @(S) real(ifft([0; Ifft_Length*frequency_spacing*S; zeros(Ifft_Length-2*no_Tones-1,1); conj(flipud(Ifft_Length*frequency_spacing*S))]));
full_synth = @(S) real(ifft([zeros(1,Nz); Ifft_Length*frequency_spacing*S; zeros(Ifft_Length-2*no_Tones-1, Nz);conj(flipud(Ifft_Length*frequency_spacing*S))]));

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
subplot(4,1,1);

plot(time_axis*1e9, recieved); grid on
xlim([0 tim_max]); xlabel('t (ns)'); ylabel('E (V/m)');
title('Recieved Signal in Time');

subplot(4,1,2);

plot(time_axis*1e9, transmitted); grid on
xlim([0 tim_max]); xlabel('t (ns)'); ylabel('E (V/m)');
title('Pulse');

subplot(4,1,3);

plot(time_axis*1e9, response_signal); grid on
xlim([0 tim_max]); xlabel('t (ns)'); ylabel('E at Rx (V/m)');
title('Recieved Conjugated with pulse');

subplot(4,1,4);

plot(time_axis*1e9, transmitted - response_signal); grid on
xlim([0 tim_max]); xlabel('t (ns)'); ylabel('E at Rx (V/m)');
title('Recieved Conjugated with pulse Minus the pulse ');


figure;
imagesc(z, time_axis*1e9, Etz);
axis xy;
colormap(jet);
colorbar;
ylim([0 tim_max]); xlabel('distance (m)'); ylabel('time (ns)');
title('E at distance/time');

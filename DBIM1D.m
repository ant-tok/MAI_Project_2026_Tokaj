clear; clc; close all;

fs_permittivity = 8.854187817e-12; 
fs_permeability = pi*4e-7;  
freq_list = (0.5:0.05:6.5)'*1e9; 
Nf = numel(freq_list);      
w = 2*pi*freq_list; 
eps_b = 40 - 9j;   
k_b = w*sqrt(fs_permeability*fs_permittivity*eps_b);  

layer_thickness = [0.002, 0.025, 0.002]; 
layer_eps = [eps_b, 54-13j, eps_b, eps_b];   
N = numel(layer_thickness) + 1;  
dist_sum = cumsum(layer_thickness);  
L = 0.04;


dz = 0.5e-3;   
Nc = round(L/dz);        
z = ((1:Nc)' - 0.5)*dz; 

eps_true = layer_eps(N)*ones(Nc, 1); 
% Masks for the layers
for n = 1:N-1   
    if n == 1
        startposition = 0; 
    else 
        startposition = dist_sum(n-1);
    end 
    eps_true((z >= startposition) & (z < dist_sum(n))) = layer_eps(n); 
end

chi_true = (eps_true - eps_b) /eps_b; 


G_D = zeros(Nc, Nc, Nf); 
G_S = zeros(Nf, Nc);  
E_i = zeros(Nc, Nf); 

for p = 1:Nf    
    s = abs(z - z.'); 

    G_D(:,:,p) = -1j*(k_b(p)*dz/2)*exp(-1j*k_b(p)*s);         
    G_D(1:Nf:end) = exp(-1j*k_b(p)*dz/2) - 1;    

    G_S(p,:) = (-1j*(k_b(p)*dz/2)*exp(-1j*k_b(p)*z)).';   
    E_i(:,p) = exp(-1j*k_b(p)*z);       
end

I = eye(Nc);

% Bassed on eqs 3 & 4
E_s = zeros(Nf, 1);       
X_fwr = diag(chi_true);             
I_fwr = eye(Nc);    
for p = 1:Nf                
    E_t = (I_fwr - G_D(:,:,p)*X_fwr) \ E_i(:,p);     
    E_s(p) = G_S(p,:)*X_fwr*E_t;
end


% Steps 1 to 7 from Magdum paper
n_iter = 14;     
lambda = 1.0e-3;
R = I;
RRE = zeros(n_iter, 1);  
MSE = zeros(n_iter, 1);     

chi = zeros(Nc, 1); 

for n = 1:n_iter    
    X = diag(chi);
    E_b = E_i;  
    G_kernel = zeros(Nf, Nc);  
    h = zeros(Nf, 1);  
    for p = 1:Nf  %freq as incidence fields
        GD_p = G_D(:,:,p);   
        GS_p = G_S(p,:);

        E_b(:,p) = (I - GD_p*X) \ E_i(:,p);

        
        G_b = GS_p / (I - X*GD_p);

        G_kernel(p,:) = G_b*diag(E_b(:,p));

        h(p) = E_s(p) - GS_p*X*E_b(:,p);  
    end

    % eqs 33 & 34 BF
    RRE(n) = norm(h)^2/norm(E_s)^2; % Relative Residual error
    MSE(n) = norm(eps_b*(1+chi) - eps_true)/norm(eps_true); 

    % BF 13
    d_chi  = (G_kernel'*G_kernel + lambda*(R'*R)) \ (G_kernel'*h);  


    chi = chi + d_chi;  

    fprintf('Iteration %2d  RRE = %.3e  MSE = %.3f\n', n, RRE(n), MSE(n));
end 


eps_rec = eps_b*(1 + chi);  % contrast to permittivity
in_mus  = (z > 0.005) & (z < 0.025); 
fprintf('Muscle layer mean eps: calced %.1f %+.1fj  true %.1f %+.1fj\n', real(mean(eps_rec(in_mus))), ...
        imag(mean(eps_rec(in_mus))), real(layer_eps(2)), imag(layer_eps(2)));

E_s_final = zeros(Nf, 1);  
X = diag(chi);
for p = 1:Nf  
    E_t = (I - G_D(:,:,p)*X) \ E_i(:,p);  
    E_s_final(p) = G_S(p,:)*X*E_t; 
end

figure;               
subplot(2,1,1);
plot(z*1e3, real(eps_true), 'k-', 'LineWidth', 1.5); hold on; grid on;
plot(z*1e3, real(eps_rec), 'r.-', 'LineWidth', 1.2);
ylabel('\epsilon'' (real)'); legend('True', 'DBIM'); title('Recovered permittivity profile');
subplot(2,1,2);
plot(z*1e3, -imag(eps_true), 'k-', 'LineWidth', 1.5); hold on; grid on;
plot(z*1e3, -imag(eps_rec), 'r.-', 'LineWidth', 1.2);
xlabel('z (mm)'); ylabel('\epsilon'''' (loss)');

figure;
semilogy(1:n_iter, RRE, 'o-', 'LineWidth', 1.2); hold on; grid on;
plot(1:n_iter, MSE, 's-', 'LineWidth', 1.2);
xlabel('Iteration'); ylabel('Error'); legend('RRE', 'MSE'); title('Convergence');
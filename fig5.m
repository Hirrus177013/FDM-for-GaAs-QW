% =========================================================================
% CHƯƠNG TRÌNH VẼ HỆ SỐ HẤP THỤ QUANG HỌC TỔNG CỘNG 
% =========================================================================
clear; clc; close all;

%% 1. THÔNG SỐ VẬT LÝ VÀ HỆ THỐNG
Nd      = 1e17;     
B_field = 0;        
T       = 300;      

L  = 15e-9;         
N  = 1500;          
dz = L / (N + 1);   
d  = 2e-9;          
z  = linspace(-L/2 + dz, L/2 - dz, N)';   

hbar = 1.0545718e-34;   
e    = 1.60217662e-19;  
eps0 = 8.85418782e-12;  
mu0  = 4 * pi * 1e-7;   
c    = 3e8;             
m0   = 9.10938356e-31;  
mstar = 0.067 * m0;      
kB   = 1.38064852e-23;  

% --- TÁCH BIỆT HAI LOẠI HẰNG SỐ ĐIỆN MÔI ---
nr       = 3.2;                  
eps_opt  = nr^2;       % Hằng số điện môi quang học (dùng cho alpha)
eps_stat = 12.58;      % Hằng số điện môi tĩnh (dùng cho Poisson)
% -------------------------------------------

tau_in  = 0.14e-12;             
gam     = hbar / tau_in;        
I_int   = 0.5 * 1e10;           

A       = (mstar * kB * T) / (pi * hbar^2);
% Tính B_const bằng eps_stat để giải Poisson
B_const = e^2 / (eps_stat * eps0); 
beta    = 1 / (kB * T);

N3d = zeros(N, 1);
N3d(abs(z) <= d/2) = Nd / d;

V0    = 228e-3 * e;    
k_nm  = 5e-9;          
V_conf = V0 * (-2 * (z/k_nm).^2 + 0.3 * (z/k_nm).^8);
vc = (2 * mstar / hbar^2) * V_conf;
b_term = (e^2 * B_field^2 / hbar^2) * z.^2;

%% 2. TẠO LƯỚI THEO TẦN SỐ GÓC OMEGA
F_array = [5e6, 15e6, 25e6]; 

E_photon_max_J = 300 * 1e-3 * e;
omega_max      = E_photon_max_J / hbar;
omega = linspace(0, omega_max, 301); 

hw_J   = hbar * omega;          
hw_meV = hw_J / e * 1000;       

alpha12_store = zeros(length(F_array), length(omega));
alpha23_store = zeros(length(F_array), length(omega));
alpha13_store = zeros(length(F_array), length(omega));

%% 3. TÍNH TOÁN HỆ SỐ HẤP THỤ alpha(omega)
Ef = 0;
VH = zeros(N, 1);

fprintf('BẮT ĐẦU TÍNH TOÁN HẤP THỤ QUANG HỌC...\n');
for step = 1:length(F_array)
    F_current = F_array(step);
    fprintf('Đang giải F = %g kV/cm... ', F_current / 1e5);
    
    f_term = (2 * mstar * e * F_current / hbar^2) * (z + L/2);
    
    for outer_iter = 1:100
        Ef_old = Ef;
        for inner_iter = 1:500
            VH_old = VH;   
            vh = (2 * mstar / hbar^2) * VH;
            
            [E_raw, psi] = schrodinger(N, b_term, f_term, vc, vh, dz);
            E_joule = E_raw * (hbar^2 / (2 * mstar));
            
            n_z = electrondensity(A, beta, Ef, E_joule, psi);
            VH_new = hatree(N, n_z, N3d, B_const, dz);
            VH = 0.5 * VH_new + (1 - 0.5) * VH_old;
            
            if max(abs(VH - VH_old))/e*1000 < 1e-3, break; end
        end
        Ef = efsolve(Nd, A, beta, E_joule, Ef_old);
        if abs(Ef - Ef_old)/e*1000 < 1e-3, break; end
    end
    
    M11 = e * sum(psi(:,1) .* z .* psi(:,1)) * dz;
    M22 = e * sum(psi(:,2) .* z .* psi(:,2)) * dz;
    M33 = e * sum(psi(:,3) .* z .* psi(:,3)) * dz;
    
    M12 = e * sum(psi(:,1) .* z .* psi(:,2)) * dz;
    M23 = e * sum(psi(:,2) .* z .* psi(:,3)) * dz;
    M13 = e * sum(psi(:,1) .* z .* psi(:,3)) * dz;
    
    % --- Dùng eps_opt cho các phương trình quang học ---
    pre1 = omega * sqrt(mu0 / (eps_opt * eps0));
    pre3 = -2 * omega * sqrt(mu0 / (eps_opt * eps0)) * (I_int / (eps0 * nr * c));
    % ---------------------------------------------------
    
    pairs = [1 2; 2 3; 1 3];
    for p = 1:3
        i = pairs(p, 1);
        f = pairs(p, 2);
        
        dE = E_joule(f) - E_joule(i);
        
        if p == 1,     M_if = M12; M_ii = M11; M_ff = M22;
        elseif p == 2, M_if = M23; M_ii = M22; M_ff = M33;
        else,          M_if = M13; M_ii = M11; M_ff = M33;
        end
        
        sigma_if = (mstar * kB * T) / (pi * hbar^2 * L) * ...
                   log( (1 + exp((Ef - E_joule(i))/(kB*T))) ./ (1 + exp((Ef - E_joule(f))/(kB*T))) );
               
        term_denom = (dE - hw_J).^2 + gam^2; 
        
        alpha1 = pre1 .* (abs(M_if)^2 * sigma_if * gam) ./ term_denom;
        
        bracket3 = 1 - (abs(M_ff - M_ii)^2 / (4 * abs(M_if)^2)) .* ...
                   ( (dE - hw_J).^2 - gam^2 + 2*dE*(dE - hw_J) ) ./ (dE^2 + gam^2);
        alpha3 = pre3 .* (abs(M_if)^4 * sigma_if * gam) ./ (term_denom.^2) .* bracket3;
        
        alpha_tot_cm = (alpha1 + alpha3) / 100;
        
        if p == 1,     alpha12_store(step, :) = alpha_tot_cm;
        elseif p == 2, alpha23_store(step, :) = alpha_tot_cm;
        else,          alpha13_store(step, :) = alpha_tot_cm;
        end
    end
    fprintf('Xong!\n');
end

%% 4. VẼ ĐỒ THỊ 
figure('Color', 'w', 'Position', [100, 50, 800, 900]);

styles = {'-', '--', ':'}; 
colors = {'b', 'k', 'r'};  
labels = {'F = 50 kV cm^{-1}', 'F = 150 kV cm^{-1}', 'F = 250 kV cm^{-1}'};

% Subplot a) alpha_12
subplot(3, 1, 1); hold on; box on; grid on;
for s = 1:3
    plot(hw_meV, alpha12_store(s, :), 'Color', colors{1}, 'LineStyle', styles{s}, 'LineWidth', 2.5, 'DisplayName', labels{s});
end
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridLineStyle', '--');
xlim([0 300]);
ylabel('Total optical absorption (\alpha_{12}) (cm^{-1})', 'FontWeight', 'bold');
legend('Location', 'northeast');
text(10, max(alpha12_store(:))*0.85, 'a) \alpha_{12}', 'FontSize', 12, 'FontWeight', 'bold');

% Subplot b) alpha_23
subplot(3, 1, 2); hold on; box on; grid on;
for s = 1:3
    plot(hw_meV, alpha23_store(s, :), 'Color', colors{2}, 'LineStyle', styles{s}, 'LineWidth', 2.5, 'DisplayName', labels{s});
end
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridLineStyle', '--');
xlim([0 300]);
ylabel('Total optical absorption (\alpha_{23}) (cm^{-1})', 'FontWeight', 'bold');
legend('Location', 'northeast');
text(10, max(alpha23_store(:))*0.85, 'b) \alpha_{23}', 'FontSize', 12, 'FontWeight', 'bold');

% Subplot c) alpha_13
subplot(3, 1, 3); hold on; box on; grid on;
for s = 1:3
    plot(hw_meV, alpha13_store(s, :), 'Color', colors{3}, 'LineStyle', styles{s}, 'LineWidth', 2.5, 'DisplayName', labels{s});
end
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridLineStyle', '--');
xlim([0 300]);
xlabel('Photon energy (meV)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Total optical absorption (\alpha_{13}) (cm^{-1})', 'FontWeight', 'bold');
legend('Location', 'northeast');
text(10, max(alpha13_store(:))*0.85, 'c) \alpha_{13}', 'FontSize', 12, 'FontWeight', 'bold');

sgtitle('Variation of \alpha_{12}, \alpha_{23}, and \alpha_{13} as a function of incident photon energy', 'FontWeight', 'bold', 'FontSize', 14);
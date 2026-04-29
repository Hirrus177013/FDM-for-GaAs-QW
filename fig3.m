% =========================================================================
% CHƯƠNG TRÌNH QUÉT ĐIỆN TRƯỜNG F VÀ VẼ ĐỘ LỆCH NĂNG LƯỢNG
% (Giải 6 mức, xuất đồ thị 4 mức thấp nhất)
% =========================================================================
clear; clc; close all;

%% 1. CÀI ĐẶT THÔNG SỐ
Nd      = 1e17;  
B_field = 0;     
Nstate  = 4;     % Chỉ lấy 4 mức để vẽ đồ thị

F_array = linspace(5e6, 25e6, 41);
E_results = zeros(length(F_array), Nstate); 

%% 2. KÍCH THƯỚC LƯỚI
L  = 15e-9;         
N  = 1500;          
dz = L / (N + 1);   
d  = 2e-9;          
z  = linspace(-L/2 + dz, L/2 - dz, N)';   

%% 3. HẰNG SỐ VẬT LÝ
hbar  = 1.0545718e-34;   
e     = 1.60217662e-19;  
eps0  = 8.85418782e-12;  
nr    = 3.2;             
eps   = 12.58;            
m0    = 9.10938356e-31;  
mstar = 0.067 * m0;      
kB    = 1.38064852e-23;  
T     = 300;             

%% 4. BIẾN TRUNG GIAN
A       = (mstar * kB * T) / (pi * hbar^2);
B_const = e^2 / (eps * eps0);
beta    = 1 / (kB * T);

%% 5. MẬT ĐỘ DONOR N3d(z)
N3d = zeros(N, 1);
N3d(abs(z) <= d/2) = Nd / d;

%% 6. THẾ NĂNG TĨNH
V0    = 228e-3 * e;    
k_nm  = 5e-9;          
beta1 = -2;
beta2 = 0.3;
V_conf = V0 * (beta1 * (z/k_nm).^2 + beta2 * (z/k_nm).^8);

vc = (2 * mstar / hbar^2) * V_conf;
b_term = (e^2 * B_field^2 / hbar^2) * z.^2;

%% 7. QUÉT ĐIỆN TRƯỜNG F
max_outer = 100;    
max_inner = 500;    
tol_Ef    = 1e-3;   
tol_VH    = 1e-3;   
lambda    = 0.5;    

fprintf('\n=====================================================\n');
fprintf('  BẮT ĐẦU QUÉT ĐIỆN TRƯỜNG F (%d ĐIỂM)\n', length(F_array));
fprintf('=====================================================\n');

Ef = 0;
VH = zeros(N, 1);

for step = 1:length(F_array)
    F_current = F_array(step);
    fprintf('Đang giải F = %g V/m (Bước %d/%d)... ', F_current, step, length(F_array));
    
    f_term = (2 * mstar * e * F_current / hbar^2) * (z + L/2);
    
    for outer_iter = 1:max_outer
        Ef_old = Ef;
        
        for inner_iter = 1:max_inner
            VH_old = VH;   
            vh = (2 * mstar / hbar^2) * VH;
            
            [E_raw, psi] = schrodinger(N, b_term, f_term, vc, vh, dz);
            E_joule = E_raw * (hbar^2 / (2 * mstar));
            
            if Nd == 0
                n = zeros(N, 1);
            else
                n = electrondensity(A, beta, Ef, E_joule, psi);
            end
            
            VH_new = hatree(N, n, N3d, B_const, dz);
            VH = lambda * VH_new + (1 - lambda) * VH_old;
            
            dVH_meV = max(abs(VH - VH_old)) / e * 1000;
            if dVH_meV < tol_VH
                break;
            end
        end
        
        Ef = efsolve(Nd, A, beta, E_joule, Ef_old);
        dEf_meV = abs(Ef - Ef_old) / e * 1000;
        
        if dEf_meV < tol_Ef
            break;
        end
    end
    
    E_results(step, :) = E_joule(1:Nstate) / e * 1000;
    fprintf('Xong! (Lặp: %d outer)\n', outer_iter);
end

%% 8. VẼ ĐỒ THỊ
F_kVcm = F_array / 1e5; 

dE_21 = E_results(:, 2) - E_results(:, 1); 
dE_31 = E_results(:, 3) - E_results(:, 1); 
dE_32 = E_results(:, 3) - E_results(:, 2); 

figure('Color', 'w', 'Position', [150, 150, 700, 500]);
hold on; box on; grid on;

set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on', 'GridLineStyle', '--');

plot(F_kVcm, dE_21, 'b-',  'LineWidth', 2.5, 'DisplayName', 'E_2 - E_1');
plot(F_kVcm, dE_31, 'r:',  'LineWidth', 2.5, 'DisplayName', 'E_3 - E_1');
plot(F_kVcm, dE_32, 'k--', 'LineWidth', 2.5, 'DisplayName', 'E_3 - E_2');

xlim([50 250]);
ylim([0 250]);
xticks(50:50:250);

xlabel('F (kV.cm^{-1})', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Energy separation (meV)', 'FontWeight', 'bold', 'FontSize', 12);

legend('Location', 'east', 'FontSize', 11);

text(60, 230, 'B = 0 T', 'FontSize', 11, 'FontWeight', 'bold');
text(60, 215, sprintf('N_d = 10^{%g} m^{-2}', log10(Nd)), 'FontSize', 11, 'FontWeight', 'bold');

hold off;
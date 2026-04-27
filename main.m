% ==========================================
% FILE: main.m
% MÔ PHỎNG GIẾNG LƯỢNG TỬ GaAs PHA TẠP
% ==========================================
clear; clc; close all;

% --- 1. KHU VỰC SETTING (DỄ DÀNG CHỈNH SỬA) ---
Nd = 1e17;        % Lượng donor trên m^2 [cite: 12]
B_mag = 0;         % Từ trường B (Tesla) [cite: 13]
F_field = 0;       % Điện trường F (V/m) [cite: 14]
Nstate = 4;        % Số mức năng lượng quan tâm sau cùng [cite: 15]

% --- 2. THÔNG SỐ HỆ THỐNG KÍCH THƯỚC ---
L = 15e-9;         % Chiều rộng cấu trúc (m) [cite: 7]
N = 1500;          % Số điểm chia lưới [cite: 8]
dz = L / (N + 1);  % Bước không gian [cite: 9]
d = 2e-9;          % Chiều rộng lớp doped (m) [cite: 10]

% --- 3. CÁC HẰNG SỐ VẬT LÝ ---
hbar = 1.0545718e-34;  % Hằng số Planck rút gọn (J.s) [cite: 18]
e = 1.60217663e-19;    % Điện tích cơ bản (C) [cite: 18]
eps0 = 8.8541878e-12;  % Độ thẩm thấu chân không (F/m) [cite: 18]
nr = 3.2;              % Chiết suất [cite: 18]
eps = nr^2;            % Hằng số điện môi tương đối [cite: 18]
m0 = 9.1093837e-31;    % Khối lượng electron (kg) [cite: 18]
mstar = 0.067 * m0;    % Khối lượng hiệu dụng e [cite: 18]
kB = 1.380649e-23;     % Hằng số Boltzmann (J/K) [cite: 18]
T = 300;               % Nhiệt độ (K) [cite: 18]

% Các biến đã được đặt [cite: 19]
A = (mstar * kB * T) / (pi * hbar^2); % [cite: 20]
B= e^2 / (eps * eps0);         % Hằng số Poisson 
beta = 1 / (kB * T);                  % [cite: 22]

% --- 4. TẠO LƯỚI KHÔNG GIAN VÀ MẬT ĐỘ DONOR ---
% Lưới z trải từ -L/2 đến L/2
z = linspace(-L/2 + dz, L/2 - dz, N)'; 

% Mật độ Donor N3d [cite: 16]
N3d = zeros(N, 1);
N3d(abs(z) <= d/2) = Nd / d; 

% --- 5. TÍNH CÁC LOẠI THẾ (POTENTIALS) CHƯA BỊ NHIỄU ---
% Thế giam hãm V_conf [cite: 24]
V0 = 228e-3 * e; % 228 meV đổi sang Joule
k_param = 5e-9;
beta1 = -2;
beta2 = 0.3;
V_conf = V0 * (beta1 * (z/k_param).^2 + beta2 * (z/k_param).^8);

% Các thành phần đưa vào phương trình Schrodinger [cite: 25, 26, 27]
vc = (2 * mstar / hbar^2) * V_conf;
b = (e^2 * B_mag^2 / hbar^2) * z.^2;
f = (2 * mstar * e * F_field / hbar^2) * (z + L/2);

% --- 6. VÒNG LẶP TỰ NHẤT QUÁN (SELF-CONSISTENCY) ---
VH = zeros(N, 1); % [cite: 31]
Ef_old = 0;       
tolerance = 1e-3; % Biến thiên Ef < 10^-3 meV [cite: 37]
lambda = 0.1;     % Hệ số trộn vH để hội tụ (phương trình 9) 

fprintf('Bắt đầu giải hệ tự nhất quán Schrodinger-Poisson...\n');
for iter = 1:500 % [cite: 32]
    % Tính vh đưa vào Schrodinger [cite: 28]
    vh = (2 * mstar / hbar^2) * VH;
    
    % Giải phương trình Schrodinger [cite: 33]
    [E_scaled, psi] = schrodinger(N, b, f, vc, vh, dz);
    
    % Trả E_scaled về lại đơn vị Joule
    E_joule = E_scaled * (hbar^2 / (2 * mstar));
    
    % Giải mức Fermi [cite: 34]
    Ef = efsolve(Nd, A, beta, E_joule);
    
    % Kiểm tra biến thiên Ef (đổi sang meV) [cite: 37]
    dEf_meV = abs(Ef - Ef_old) / e * 1000;
    if iter > 1 && dEf_meV < tolerance
        fprintf('-> Hội tụ thành công sau %d vòng lặp! (dEf = %.2e meV)\n', iter, dEf_meV);
        break;
    end
    Ef_old = Ef;
    
    % Giải mật độ n(z) [cite: 35]
    n = electrondensity(A, beta, Ef, E_joule, psi);
    
    % Giải lại thế Hartree VH(z) [cite: 36]
    % Gọi hàm hatree với B_const theo cập nhật mới nhất của chúng ta
    VH_new = hatree(N, n, N3d, B, dz);
    
    % Trộn nghiệm VH để tránh dao động phân kỳ 
    VH = lambda * VH_new + (1 - lambda) * VH;
    
    if mod(iter, 10) == 0
        fprintf('  Vòng %d: dEf = %.2e meV\n', iter, dEf_meV);
    end
end

% --- 7. OUTPUT KẾT QUẢ VÀ ĐỒ THỊ ---
% Đổi năng lượng sang meV và z sang nm
Ef_meV = Ef / e * 1000;
E_meV = E_joule / e * 1000;
% --- Phục hồi các thế năng từ b và f ---
V_mag = b * (hbar^2 / (2 * mstar));    % Thế năng do từ trường B
V_elec = f * (hbar^2 / (2 * mstar));   % Thế năng do điện trường F

% --- Tổng hợp TẤT CẢ các thế năng ---
V_tot_physical = V_conf + VH + V_mag + V_elec;
V_tot_meV = V_tot_physical / e * 1000;
z_nm = z * 1e9;

% In ra Command Window [cite: 39, 40, 41]
fprintf('\n=================================\n');
fprintf('Mức Fermi Ef = %.4f meV\n', Ef_meV);
fprintf('%d mức năng lượng thấp nhất:\n', Nstate);
for i = 1:Nstate
    fprintf('  E_%d = %.4f meV\n', i, E_meV(i));
end
fprintf('=================================\n');

% Vẽ đồ thị theo đúng style tham khảo [cite: 42, 47]
figure('Color', 'w', 'Position', [100, 100, 750, 500]);
hold on; box on; grid on;

% Cài đặt mảng màu và đường vẽ cho các mức
colors = {'r', 'g', 'k', 'm'};
lines  = {'-', '--', '--', '-.'};
names  = {'Ground state', '1^{st} excited', '2^{nd} excited', '3^{rd} excited'};

% Vẽ các hàm sóng được đẩy lên mức năng lượng E_i
scale_factor = 300 / max(max(psi(:, 1:Nstate).^2));

for i = 1:Nstate
    psi2_plot = E_meV(i) + scale_factor * (psi(:, i).^2);
    plot(z_nm, psi2_plot, 'Color', colors{i}, 'LineStyle', lines{i}, ...
         'LineWidth', 2.5, 'DisplayName', names{i});
end

% Vẽ mức Fermi (vàng, chấm)
plot(z_nm, Ef_meV * ones(N, 1), 'Color', [0.9290 0.6940 0.1250], ...
     'LineStyle', ':', 'LineWidth', 3, 'DisplayName', 'Fermi level');

% Vẽ thế năng giam hãm (xanh đậm, liền)
plot(z_nm, V_tot_meV, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Confining potential');

% Định dạng text, trục và giới hạn
xlabel('z (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy(meV)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-7.5 7.5]);
ylim([-500 1000]);
set(gca, 'FontSize', 11, 'LineWidth', 1.2);

% Thêm chú thích text góc trái [cite: 47]
if Nd == 0
    Nd_str = 'N_d = 0 m^{-2}';
else
    Nd_str = sprintf('N_d = 10^{%d} m^{-2}', round(log10(Nd)));
end

% Các dòng text in lên đồ thị giữ nguyên
text(-6.5, 850, 'b)', 'FontSize', 14, 'FontWeight', 'bold');
text(-5.5, 850, Nd_str, 'FontSize', 13, 'FontWeight', 'bold');
text(-6.5, 850, 'b)', 'FontSize', 14, 'FontWeight', 'bold');
text(-5.5, 850, Nd_str, 'FontSize', 13, 'FontWeight', 'bold');
% --- Xử lý định dạng chữ cho Điện trường F đẹp mắt ---
if F_field == 0
    F_str = '0';
else
    % Tách phần hệ số và phần mũ
    F_exp = floor(log10(abs(F_field)));
    F_base = F_field / (10^F_exp);
    
    % Nếu hệ số là 1 (vd: 10^6), chỉ in 10^6. Nếu khác 1 (vd: 8*10^6), in 8 \times 10^6
    if F_base == 1
        F_str = sprintf('10^{%d}', F_exp);
    else
        F_str = sprintf('%g \\times 10^{%d}', F_base, F_exp);
    end
end

% Ghép vào chuỗi text cuối cùng
BF_str = sprintf('B = %g T; F = %s V/m;', B_mag, F_str);

% In lên đồ thị (vẫn giữ nguyên tọa độ)
text(-5.5, 750, BF_str, 'FontSize', 13, 'FontWeight', 'bold');
% Cài đặt Legend
legend('Location', 'northeast', 'FontSize', 11);
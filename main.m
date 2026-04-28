% =========================================================================
% CHƯƠNG TRÌNH CHÍNH: MAIN.M
% Mô phỏng mức năng lượng trong giếng lượng tử GaAs pha tạp
% =========================================================================
clear; clc; close all;

%% 1. CÀI ĐẶT 4 THÔNG SỐ HỆ THỐNG  <-- Chỉnh tại đây
Nd     = 1e17;  % Lượng donor trên m^2      (đặt 0 nếu không pha tạp)
B_field = 0;    % Từ trường B (Tesla)
F_field = 0;    % Điện trường F (V/m)
Nstate  = 4;    % Số mức năng lượng quan tâm sau cùng

%% 2. CÁC THÔNG SỐ HỆ (KÍCH THƯỚC)
L  = 15e-9;         % Chiều rộng giếng thế (m)
N  = 1500;          % Số điểm lưới (giảm xuống 500 nếu máy chậm)
dz = L / (N + 1);   % Khoảng cách giữa các điểm lưới
d  = 2e-9;          % Chiều rộng lớp doped (m)

z = linspace(-L/2 + dz, L/2 - dz, N)';   % Vector vị trí (cột)

%% 3. CÁC HẰNG SỐ VẬT LÝ
hbar = 1.0545718e-34;   % J.s
e    = 1.60217662e-19;  % C
eps0 = 8.85418782e-12;  % F/m
nr   = 3.2;             % Chiết suất GaAs
eps  = nr^2;            % Hằng số điện môi tương đối
m0   = 9.10938356e-31;  % kg
mstar = 0.067 * m0;     % Khối lượng hiệu dụng electron GaAs
kB   = 1.38064852e-23;  % J/K
T    = 300;             % Nhiệt độ (K)

%% 4. CÁC BIẾN ĐÃ ĐẶT: A, B_const, beta
% A = (mstar * kB * T) / (pi * hbar^2)
A      = (mstar * kB * T) / (pi * hbar^2);
% B_const = e^2 / (eps * eps0)  [dùng trong hatree.m]
B_const = e^2 / (eps * eps0);
% beta = 1 / (kB * T)
beta   = 1 / (kB * T);

%% 5. MẬT ĐỘ DONOR N3d(z)
% N3d = Nd/d khi |z| <= d/2, bằng 0 bên ngoài
N3d = zeros(N, 1);
N3d(abs(z) <= d/2) = Nd / d;

%% 6. CÁC LOẠI THẾ NĂNG TĨNH (đơn vị được chuẩn hoá cho schrodinger.m)
% --- Thế giam hãm dị hướng (anharmonic confinement) ---
V0    = 228e-3 * e;    % Biên độ thế (Joule)
k_nm  = 5e-9;          % Tham số tỉ lệ (m)
beta1 = -2;
beta2 = 0.3;
V_conf = V0 * (beta1 * (z/k_nm).^2 + beta2 * (z/k_nm).^8);

% --- Chuyển sang dạng chuẩn hoá cho Schrodinger solver ---
% vc = (2*mstar/hbar^2) * V_conf
vc = (2 * mstar / hbar^2) * V_conf;

% --- Thành phần từ trường: b = e^2*B^2/hbar^2 * z^2 ---
b_term = (e^2 * B_field^2 / hbar^2) * z.^2;

% --- Thành phần điện trường: f = -2*mstar*e*F/hbar^2 * (z + L/2) ---
f_term = (2 * mstar * e * F_field / hbar^2) * (z + L/2);

%% 7. VÒNG LẶP TỰ NHẤT QUÁN (SELF-CONSISTENT LOOP)
% Tham số hội tụ
max_outer = 100;    % Số lần lặp tối đa (outer)
max_inner = 500;    % Số lần lặp tối đa (inner)
tol_Ef    = 1e-3;   % Ngưỡng hội tụ Ef (meV)
tol_VH    = 1e-3;   % Ngưỡng hội tụ VH (meV)
lambda    = 0.5;   % Hệ số mixing VH (giảm dao động)

fprintf('\n=====================================================\n');
fprintf('  BẮT ĐẦU VÒNG LẶP TỰ NHẤT QUÁN\n');
fprintf('=====================================================\n');

% --- Khởi tạo Ef = 0 và VH = 0 (một lần duy nhất trước outer loop) ---
Ef = 0;
VH = zeros(N, 1);

for outer_iter = 1:max_outer
    Ef_old = Ef;

    % ===== INNER LOOP: hội tụ VH (VH được kế thừa từ vòng outer trước) =====
    for inner_iter = 1:max_inner
        VH_old = VH;   % Lưu VH hiện tại làm reference cho mixing

        % Bước 1: Chuyển VH sang dạng chuẩn hoá
        vh = (2 * mstar / hbar^2) * VH;

        % Bước 2: Giải phương trình Schrödinger
        % [E, psi] = schrodinger(N, b, f, vc, vh, dz)
        [E_raw, psi] = schrodinger(N, b_term, f_term, vc, vh, dz);

        % Bước 3: Phục hồi đơn vị năng lượng (Joule)
        E_joule = E_raw * (hbar^2 / (2 * mstar));

        % Bước 4: Tính mật độ electron n(z)
        if Nd == 0
            n = zeros(N, 1);
        else
            n = electrondensity(A, beta, Ef, E_joule, psi);
        end

        % Bước 5: Giải Poisson để lấy VH mới, rồi mixing
        % hatree(N, n, N3d, B_const, dz)  [5 tham số theo docx]
        VH_new = hatree(N, n, N3d, B_const, dz);
        VH = lambda * VH_new + (1 - lambda) * VH_old;

        % Kiểm tra hội tụ inner (đơn vị meV)
        dVH_meV = max(abs(VH - VH_old)) / e * 1000;
        if dVH_meV < tol_VH
            break;
        end
    end

    % ===== OUTER LOOP: cập nhật Ef =====
    % efsolve(Nd, A, beta, E, Ef_old)  [5 tham số]
    Ef = efsolve(Nd, A, beta, E_joule, Ef_old);

    dEf_meV = abs(Ef - Ef_old) / e * 1000;
    fprintf('[Outer %3d]  Ef = %8.3f meV  |  dEf = %.4e meV  |  Inner steps: %d\n', ...
            outer_iter, Ef/e*1000, dEf_meV, inner_iter);

    % Điều kiện dừng toàn cục: chỉ dừng khi dEf thực sự nhỏ hơn tolerance
    if dEf_meV < tol_Ef
        fprintf('\n==> HỘI TỤ THÀNH CÔNG sau %d vòng lặp ngoài!\n', outer_iter);
        break;
    end
    if outer_iter == max_outer
        fprintf('\n==> CẢNH BÁO: Chưa hội tụ sau %d vòng lặp!\n', max_outer);
    end
end

%% 8. OUTPUT KẾT QUẢ
% --- Tính thế năng tổng cộng ---
V_mag      = b_term  * (hbar^2 / (2 * mstar));   % Thành phần từ trường (J)
V_elec     = f_term  * (hbar^2 / (2 * mstar));   % Thành phần điện trường (J)
V_tot_meV  = (V_conf + VH + V_mag + V_elec) / e * 1000;  % meV

% --- Các mức năng lượng ---
E_meV  = E_joule / e * 1000;   % meV
Ef_meV = Ef      / e * 1000;   % meV

fprintf('\nKẾT QUẢ CUỐI CÙNG:\n');
fprintf('  Mức Fermi Ef = %.4f meV\n', Ef_meV);
fprintf('  %d mức năng lượng thấp nhất:\n', Nstate);
for i = 1:Nstate
    fprintf('  E%d = %.4f meV\n', i, E_meV(i));
end

%% 9. VẼ ĐỒ THỊ
figure('Color', 'w', 'Position', [100, 100, 800, 500]);
hold on; box on; grid on;

z_nm = z * 1e9;  % Đổi sang nm

% Tỉ lệ vẽ hàm sóng (|psi|^2 được scale để hiển thị rõ)
scale = 300 / max(max(psi(:, 1:Nstate).^2));

colors     = {'r',  'g',   'k',    'm'};
linestyles = {'-',  '--',  '--',   '-.'};
labels     = {'Ground state', '1^{st} excited', '2^{nd} excited', '3^{rd} excited'};

% Vẽ |psi|^2 dịch chuyển lên mức năng lượng tương ứng
for i = 1:Nstate
    psi2_shifted = E_meV(i) + scale * psi(:, i).^2;
    plot(z_nm, psi2_shifted, ...
         'Color', colors{i}, 'LineStyle', linestyles{i}, ...
         'LineWidth', 2.5, 'DisplayName', labels{i});
end

% Vẽ mức Fermi
yline(Ef_meV, 'Color', [0.9290 0.6940 0.1250], ...
      'LineStyle', ':', 'LineWidth', 2.5, 'DisplayName', 'Fermi level');

% Vẽ thế năng tổng cộng
plot(z_nm, V_tot_meV, 'b-', 'LineWidth', 2.5, 'DisplayName', 'V_{total}');

% Thiết lập trục và nhãn
xlabel('z (nm)', 'FontWeight', 'bold', 'FontSize', 12);
ylabel('Energy (meV)', 'FontWeight', 'bold', 'FontSize', 12);
xlim([-7.5 7.5]);
ylim([-500 1000]);
legend('Location', 'northeast', 'FontSize', 10);

% --- Hiển thị thông số Nd, B, F lên hình ---
% Định dạng Nd
if Nd == 0
    Nd_str = 'N_d = 0 m^{-2}';
else
    Nd_str = sprintf('N_d = 10^{%g} m^{-2}', log10(Nd));
end

% Định dạng F
if F_field == 0
    F_str = '0';
else
    F_exp  = floor(log10(abs(F_field)));
    F_base = F_field / 10^F_exp;
    F_str  = sprintf('%g \\times 10^{%d}', F_base, F_exp);
end

text(-7, 850, sprintf('b)     %s', Nd_str), ...
     'FontSize', 12, 'FontWeight', 'bold');
text(-7, 700, sprintf('         B = %g T;  F = %s V/m', B_field, F_str), ...
     'FontSize', 12, 'FontWeight', 'bold');

hold off;
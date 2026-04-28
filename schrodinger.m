function [E, psi] = schrodinger(N, b, f, vc, vh, dz)
    % Hàm giải phương trình Schrodinger 1D bằng sai phân hữu hạn
    % TỐI ƯU CỰC HẠN: Giải 6 mức + Chuẩn hóa Grid + Kỹ thuật Shift

    if nargin < 6, dz = 1; end

    % 1. Đảm bảo các thế năng là vector cột
    b = b(:); f = f(:); vc = vc(:); vh = vh(:);
    V_total = b + f + vc + vh;

    % 2. Xây dựng Hamiltonian chuẩn hóa (nhân dz^2) để triệt tiêu số siêu lớn
    e_vec = ones(N, 1);
    T_grid = spdiags([-e_vec, 2*e_vec, -e_vec], [-1, 0, 1], N, N); % Động năng
    V_grid = spdiags(V_total .* (dz^2), 0, N, N);                 % Thế năng
    H_grid = T_grid + V_grid;
    
    % Ép đối xứng tuyệt đối
    H_grid = (H_grid + H_grid') / 2;

    % 3. Kỹ thuật Shift & eigs cho 6 mức năng lượng
    num_states = 6;  % CHỈ LẤY ĐÚNG 6 MỨC
    opts.disp = 0;   % Tắt thông báo rác
    
    % Tìm sigma (điểm trũng nhất của giếng thế) làm mốc neo cho hàm eigs
    sigma = min(V_total .* (dz^2)); 
    [V_eig, D_eig] = eigs(H_grid, num_states, sigma, opts);

    % 4. Phục hồi đơn vị gốc
    E_grid = diag(D_eig);
    [E_grid, sort_idx] = sort(E_grid);
    E = E_grid / (dz^2);
    
    psi = V_eig(:, sort_idx);
    psi = psi / sqrt(dz);
end
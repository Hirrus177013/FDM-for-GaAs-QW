function [E, psi] = schrodinger(N, b, f, vc, vh, dz)
    % Kiểm tra đầu vào để tìm thủ phạm gây NaN/Inf
    vars = {b, f, vc, vh, dz};
    names = {'b (Từ trường)', 'f (Điện trường)', 'vc (Giam hãm)', 'vh (Hartree)', 'dz (Lưới)'};
    for i = 1:5
        if any(isnan(vars{i}(:))) || any(isinf(vars{i}(:)))
            error('LỖI: Biến %s chứa giá trị NaN hoặc Inf!', names{i});
        end
    end

    % Đảm bảo các thế năng là vector cột
    V_total = b(:) + f(:) + vc(:) + vh(:);

    % Xây dựng Hamiltonian chuẩn hóa (nhân dz^2) để các số hạng về cỡ O(1)
    e_vec = ones(N, 1);
    T_grid = spdiags([-e_vec, 2*e_vec, -e_vec], [-1, 0, 1], N, N); % Động năng
    V_grid = spdiags(V_total .* (dz^2), 0, N, N);                 % Thế năng
    H_grid = T_grid + V_grid;
    
    % Ép đối xứng để eigs chạy ổn định hơn
    H_grid = (H_grid + H_grid') / 2;

    % Giải 15 mức thấp nhất
    num_states = 15;
    opts.disp = 0;
    
    % Kỹ thuật SHIFT: Tìm các trị riêng quanh mức thế năng thấp nhất để tránh "badly conditioned"
    sigma = min(V_total .* (dz^2)); 
    [V_eig, D_eig] = eigs(H_grid, num_states, sigma, opts);

    % Chuyển đổi ngược lại đơn vị Joule
    E_grid = diag(D_eig);
    [E_grid, sort_idx] = sort(E_grid);
    E = E_grid / (dz^2);
    psi = V_eig(:, sort_idx);

    % Chuẩn hóa hàm sóng
    psi = psi / sqrt(dz);
end
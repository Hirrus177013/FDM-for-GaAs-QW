function n = electrondensity(A, beta, Ef, E, psi)
    % electrondensity: Tính phân bố mật độ electron n(z) trong giếng lượng tử
    %
    % ĐẦU VÀO:
    %   A    : Hệ số mật độ trạng thái (scalar)
    %   beta : 1 / (kB * T) (scalar)
    %   Ef   : Mức năng lượng Fermi (scalar)
    %   E    : Vector các mức năng lượng E_i (kích thước M x 1)
    %   psi  : Ma trận hàm sóng (kích thước N x M), 
    %          trong đó N là số điểm z, M là số mức năng lượng.
    % ĐẦU RA:
    %   n    : Vector phân bố mật độ electron theo không gian z (kích thước N x 1)

    % Đảm bảo E là vector cột
    E = E(:);

    % =========================================================
    % NHIỆM VỤ 1: Tính vector ntot (Tổng số electron tại mỗi mức i)
    % CÔNG THỨC GỐC: ntot = A * log(1 + exp(beta * (Ef - E)))
    % =========================================================
    
    u = beta * (Ef - E);
    
    % Sử dụng Stable Softplus trick để tránh lỗi tràn số (Overflow)
    ntot = A * (max(u, 0) + log(1 + exp(-abs(u))));
    
    % Đảm bảo ntot là vector cột (kích thước M x 1)
    ntot = ntot(:);

    % =========================================================
    % NHIỆM VỤ 2: Tính phân bố mật độ n(z)
    % =========================================================

    n = (psi.^2) * ntot;
    
end
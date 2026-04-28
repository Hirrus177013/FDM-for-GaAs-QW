function n = electrondensity(A, beta, Ef, E, psi)
    % electrondensity: Tính phân bố mật độ electron n(z)
    % Đã tự động tương thích với đầu vào E (15x1) và psi (1500x15)

    % Đảm bảo E là vector cột
    E = E(:);

    % =========================================================
    % NHIỆM VỤ 1: Tính vector ntot (Tổng số electron tại mỗi mức i)
    % =========================================================
    u = beta * (Ef - E);
    
    % Sử dụng Stable Softplus trick để tránh lỗi tràn số (Overflow)
    ntot = A * (max(u, 0) + log(1 + exp(-abs(u))));
    
    % Đảm bảo ntot là vector cột (kích thước 15 x 1)
    ntot = ntot(:);

    % =========================================================
    % NHIỆM VỤ 2: Tính phân bố mật độ n(z)
    % =========================================================
    % psi.^2 là (1500x15), ntot là (15x1). Kết quả n là (1500x1)
    n = (psi.^2) * ntot;
end
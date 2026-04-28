function n = electrondensity(A, beta, Ef, E, psi)
    % Tính phân bố mật độ electron n(z)
    % Tự động tương thích với ma trận 6 mức

    E = E(:);
    u = beta * (Ef - E);
    
    % Sử dụng Stable Softplus trick chống tràn số
    ntot = A * (max(u, 0) + log(1 + exp(-abs(u))));
    ntot = ntot(:);

    % Phép nhân ma trận tốc độ cao: (1500x6) * (6x1) = (1500x1)
    n = (psi.^2) * ntot;
end
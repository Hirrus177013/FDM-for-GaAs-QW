function VH = hatree(N, n, N3d, B, dz)
    % Hàm tính thế Hartree bằng thuật toán Thomas
    % Phương trình: VH''(z) = B * (N3d - n(z))
    %
    % ĐẦU VÀO:
    %   N   : Kích thước ma trận (số điểm chia lưới phần lõi)
    %   n   : Vector mật độ electron (kích thước N x 1)
    %   N3d : Mật độ khối donor (scalar hoặc vector kích thước N x 1)
    %   B   : Hằng số tỉ lệ
    %   dz  : Bước chia không gian
    % ĐẦU RA:
    %   VH  : Vector thế Hartree (kích thước N x 1)
    
    % Ép n về vector cột
    n = n(:); 
    
    % Tính vector vế phải (RHS) với biến N3d
    rho = (dz^2) * B * (N3d - n); 
    
    % Khởi tạo mảng alpha và gamma đúng kích thước N
    alpha = zeros(N, 1);
    gamma = zeros(N, 1);
    
    % Các hằng số của ma trận sai phân D
    d_j = -2; % Đường chéo chính
    a_j = 1;  % Đường chéo dưới
    c_j = 1;  % Đường chéo trên
    
    % --- BƯỚC 1: QUÉT THUẬN (Forward Sweep) ---
    % Khởi tạo tại j = 1 (biên 0)
    alpha(1) = -c_j / d_j;
    gamma(1) = rho(1) / d_j;
    
    % Quét từ j = 2 đến N
    for j = 2:N
        mau_so = d_j + a_j * alpha(j-1);
        alpha(j) = -c_j / mau_so;
        gamma(j) = (rho(j) - a_j * gamma(j-1)) / mau_so;
    end
    
    % --- BƯỚC 2: QUÉT NGƯỢC (Backward Sweep) ---
    VH = zeros(N, 1);
    
    % Khởi tạo tại j = N (biên N+1 = 0)
    VH(N) = gamma(N);
    
    % Quét lùi từ N-1 về 1
    for j = N-1:-1:1
        VH(j) = alpha(j) * VH(j+1) + gamma(j);
    end
end
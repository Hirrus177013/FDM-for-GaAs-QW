function [E, psi] = schrodinger(N, b, f, vc, vh, dz)
    % Hàm giải phương trình Schrodinger 1D bằng sai phân hữu hạn
    % Đầu ra:
    % E   : Vector cột chứa TẤT CẢ N mức năng lượng
    % psi : Ma trận chứa TẤT CẢ N hàm sóng

    if nargin < 6
        dz = 1; 
    end

    % Ép TẤT CẢ các thành phần thế năng về vector cột (đề phòng lỗi kích thước)
    b = b(:);
    f = f(:);
    vc = vc(:); 
    vh = vh(:);

    % --- ĐÃ SỬA LỖI TẠI ĐÂY ---
    % Cộng trực tiếp các vector thế năng lại với nhau
    V_total = b + f + vc + vh;

    % Xây dựng Hamiltonian
    e = ones(N, 1);
    D2 = spdiags([e, -2*e, e], [-1, 0, 1], N, N) / (dz^2);
    V_matrix = spdiags(V_total, 0, N, N);
    
    % H = -D2 + V
    H = -D2 + V_matrix;

    % Giải toàn bộ phổ bằng eig (chuyển thưa thành đặc)
    [V_eig, D_eig] = eig(full(H));

    % Trích xuất E và sắp xếp từ thấp lên cao
    E = diag(D_eig);
    [E, sort_idx] = sort(E);
    
    % Sắp xếp lại các cột hàm sóng tương ứng với E
    psi = V_eig(:, sort_idx);

    % Chuẩn hóa hàm sóng
    psi = psi / sqrt(dz);
end
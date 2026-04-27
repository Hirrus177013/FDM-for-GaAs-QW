function Ef = efsolve(Nd, A, beta, E)
    % efsolve: Giải phương trình cân bằng mật độ tìm Ef bằng Newton-Raphson tối ưu
    % Phương trình: Nd = A * sigma_i ln(1 + exp(beta * (Ef - E_i)))
    %
    % Đầu vào:
    %   Nd   : Mật độ donor (scalar)
    %   A    : Hệ số tỉ lệ, thường là (m* kB T) / (pi hbar^2) (scalar)
    %   beta : 1 / (kB*T) (scalar)
    %   E    : Vector chứa các mức năng lượng E_i
    % Đầu ra:
    %   Ef   : Mức Fermi (scalar)

    % --- 1. Tiền xử lý dữ liệu ---
    E = E(:); % Đảm bảo E là vector cột
    
    % --- 2. Dự đoán ban đầu ---
    Ef = 0; 
    
    % --- 3. Cài đặt thông số Newton-Raphson ---
    tol = 1e-12;      
    max_iter = 100;   
    
    % --- 4. Vòng lặp ---
    for i = 1:max_iter
        % Đối số u_i
        u = beta * (Ef - E);
        
        % Tính f(Ef) dùng Stable Softplus Trick
        term_f = max(u, 0) + log(1 + exp(-abs(u)));
        val_f = A * sum(term_f) - Nd;  % Nhân A và trừ Nd
        
        % Tính đạo hàm f'(Ef)
        term_df = 1 ./ (1 + exp(beta * (E - Ef)));
        val_df = A * beta * sum(term_df); % Nhân A vào đạo hàm
        
        if abs(val_df) < 1e-20
            error('Đạo hàm quá nhỏ, Newton-Raphson thất bại.');
        end
        
        % Cập nhật Ef
        delta_Ef = val_f / val_df;
        Ef = Ef - delta_Ef;
        
        % Kiểm tra hội tụ
        if abs(delta_Ef) < tol
            return; 
        end
    end
    
    warning('efsolve không hội tụ được đến độ chính xác tolerance yêu cầu trong %d vòng lặp.', max_iter);
end
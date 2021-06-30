% Time symbol-wise IFFT and CP

function [time_grid] = ifft_sym_wise(pre_grid, para)

    %Initialising frame to zero
    frame_ifft = zeros(para.N_point,para.time_sym);

    % Taking IFFT per time symbol
    for ii=1:para.time_sym
        frame_ifft(:,ii) = ifft(pre_grid(:,ii),para.N_point);
    end

    % Adding cyclic prefix
    cp_rows = para.cp_length;
    time_grid = [frame_ifft((para.N_point-cp_rows+1):para.N_point,:);frame_ifft];

end
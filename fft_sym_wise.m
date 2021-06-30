% Time symbol-wise FFT 

function [frame_dec] = fft_sym_wise(Rx_frame, para)

    % Removing CP
    frame_CP_rem = Rx_frame(para.cp_length+1:size(Rx_frame,1),:);


    % Taking FFT per time symbol
    for ii=1:para.time_sym
        frame_fft(:,ii) = fft(frame_CP_rem(:,ii),para.N_point);
    end
    
    frame_dec = frame_fft(1:para.sc_per_PRB*para.total_PRBs,:);
    
end

% Rearrangement in grid form

function [pre_grid] = form_grid(per_user_pream_seq,para,subband_number)

    pre_grid = zeros(para.total_sc, para.time_sym);
    
    % Mini grid
    mini_grid = zeros(para.num_preamble_sc,para.time_sym);
    mini_grid(1:length(per_user_pream_seq)) = per_user_pream_seq;
    
    % Full grid arrangement
    grid_start = (subband_number-1) * para.num_preamble_sc + 1; 
    grid_end = (subband_number*para.num_preamble_sc);
    pre_grid(grid_start:grid_end,:) = mini_grid;  
    
end
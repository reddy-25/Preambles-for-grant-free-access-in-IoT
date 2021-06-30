% PREAMBLE ONLY CODE, NO DATA

close all;
clear all;

%% Parameters

rng(2020);

format long;
no_iter = 10;         % No of iterations
para.users = 2;         % No of users

% Grid dimensions
para.sc_spacing = 3.75*(10^3);
para.total_PRBs = 6;        % Number of PRBs
para.sc_per_PRB = 48;       % Assuming for 3.75kHz spacing  
para.total_sc = para.sc_per_PRB*para.total_PRBs; 
para.time_sym = 7;          
para.N_point = 512;         % Assuming that there are 48 sc per PRB with a spacing of 3.75kHz
para.cp_length = 128;   

% Preamble related parameters
prob_of_succ = 0.9;
para.total_pre_req = ceil((para.users - 1)/(1-(prob_of_succ^(2/para.users ))));
paravals_1.total_pre_req = max(para.total_pre_req,1);
para.num_preamble_PRB = 6; % Set to a factor of total PRBs
para.num_choices = para.total_PRBs/para.num_preamble_PRB; %Number of sub band options for the users to pick
para.no_of_pre = ceil(para.total_pre_req/para.num_choices); 
para.num_preamble_sc = para.sc_per_PRB*para.num_preamble_PRB;


para.long_preamble = 0;
para.short_preamble = 1;
para.spread_sequence = 0;
para.auto_correlation =1;
if(para.spread_sequence==1)
    para.no_of_pre=12;
end
% Length of the preamble
%para.pre_length = 288;
% For 6 subband choices - 48 length preamble - repeated 54 times
% For 1 subband choice - 288 length preamble - repeated 9 times
para.short_pre_length = para.num_preamble_sc;

while ~isprime(para.short_pre_length)
        para.short_pre_length = para.short_pre_length - 1;
end


%para.pre_length = (para.N_point+para.cp_length)*4;
para.pre_length = para.num_preamble_sc * para.time_sym ;  

% The actual preamble length will be the nearest prime number <= this value
    while ~isprime(para.pre_length)
        para.pre_length = para.pre_length - 1;
    end

% Simulation parameters
SNR_start_dB = -20;        % in dB
SNR_end_dB = 10;           % in dB
SNR_step_dB = 5;           % in dB

% Preambles populated in time domain or freq domain
para.pre_domain = 'FD';     % Set FD for freq domain ;; Set TD for time domain
if (para.pre_domain == "TD")
    para.num_preamble_PRB = para.total_PRBs; 
    para.num_choices = 1;
end

% Include noise
add_noise = 1;              % If 1, AWGN is considered
signal_power = 1; 
noise_mean = 0;

% Include channel
include_channel = 0;        % If 1, fading channel is considered
n_tap = 1;                  % Number of taps for the fading channel

include_lteFC = 1;          % If 1, lte fading channel is considered

% LTE Fading channel
% Channel model configuration structure
chcfg.DelayProfile = 'EPA';
chcfg.NRxAnts = 1;
chcfg.DopplerFreq = 1;
chcfg.MIMOCorrelation = 'Low';
chcfg.Seed = 1;
chcfg.InitPhase = 'Random';
chcfg.ModelType = 'GMEDS';
chcfg.NTerms = 16;
chcfg.NormalizeTxAnts = 'On';
chcfg.NormalizePathGains = 'On';
sampling_rate = para.sc_spacing*para.N_point;

% Delay tolerance
% Tolerance for the delay than can be accepted. + or - these many samples
delay_tolerance = 7;

% Frequency offset range
fo_range = 200;

% Threshold range for detection purpose
% Set threshold value to 3 to get 0.1% false alarm in AWGN for FD
% Set threshold value to 8 to get 0.1% false alarm in AWGN for TD
% For checking the performance over different thresh values, populate the array
%para.threshold = [80];  
para.threshold =[5:10];
%para.threshold =[100];
%para.threshold =[3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.8, 4, 4.2,4.4, 4.5];
%para.threshold=[1500];

%% Generating premables


[preambles_set,para] = gen_preambs_common(para);
%preambles_set_fft = fft(preambles_set);


%% Simulation Loop 
  

SNR_range = SNR_start_dB:SNR_step_dB:SNR_end_dB;
store_delay_diff = [];
%delay_list_table = [];


% Initializing other parameters to zero
false_pos = zeros(length(SNR_range),para.num_choices,length(para.threshold));
coll = zeros(length(SNR_range),para.num_choices,length(para.threshold));
count_succ = zeros(length(SNR_range),length(para.threshold));

tic
for count = 1:1:no_iter
    
    SNR_count = 1;
    
    for user_count = 1:1:para.users
        
        % User randomly picking premable (and slot number)
        subband_number = randi(para.num_choices);  % User picking subband number
        subband_select(user_count) = subband_number;
        
        pre_number = user_count;%randi(para.no_of_pre);
        pream_select(user_count) = pre_number;
        per_user_pream_seq = preambles_set(:,pre_number);  % Each user's premable sequence
        para.len_per_usr_preseq = length(per_user_pream_seq);
        
        time_grid_size = [para.N_point+para.cp_length para.time_sym];
        tx_grid = zeros(time_grid_size);   % Initialising the grid
        % Rearrangement in grid form
        if (para.pre_domain == "FD")
            [pre_grid] = form_grid(per_user_pream_seq,para,subband_number);
           
            % Taking IFFT and adding CP
            time_grid = ifft_sym_wise(pre_grid, para);
            tx_grid = time_grid;
        else
            % Assuming preamble occupies the first symbols in TD
            tx_grid(1:length(per_user_pream_seq)) = per_user_pream_seq;
        end
        
        % Rearranging grid as transmission vector
        tx_signal = tx_grid(:);
        
        
        %% CHANNEL
        
        
        if (include_channel == 1)
            if(include_lteFC == 1)
                % Transmitted waveform configuration structure
                chcfg.InitTime = count/1000;    % 1ms per iteration count
                % Number of transmit antennas and waveform sampling rate
                P = size(tx_signal,2);
                chcfg.SamplingRate = para.sc_spacing*para.N_point;
                
                % Send waveform through the channel
                [tx_signal_conv,info] = lteFadingChannel(chcfg,[tx_signal; zeros(25,P)]);
                
                chan_delay = info.ChannelFilterDelay; % Storing the delay added by the channel
            else
                h = (randn(1,n_tap)+1i*randn(1,n_tap))/sqrt(2*n_tap);    % Generating channel impulse response
                % Generate again if h > -10
                while (10*log10((abs(h)).^2) <= -10)
                    h = (randn(1,n_tap)+1i*randn(1,n_tap))/sqrt(2*n_tap);
                end
                %h_value(count) = 10*log10((abs(h))^2);
                tx_signal_conv = conv(tx_signal,h);
                chan_delay = 0;
            end
        else
            tx_signal_conv = tx_signal;
            chan_delay = 0;
        end
         
        tx_signal_conv = tx_signal_conv.';
        
        % Introducing frequency offset
        freq_offset = 0;% randi([-fo_range fo_range],1,1);
        tx_signal_index = (0:length(tx_signal_conv)-1);
        alpha = (2*pi*freq_offset)/sampling_rate;
        
        tx_signal_FO = tx_signal_conv.*(exp(1i*alpha*tx_signal_index));
        
        % Introducing time delay
        % User picks the delay randomly within the CP length
        delay_assigned(user_count) = 0; % randi([0 para.cp_length], 1, 1);
        delay_zeroes = zeros(1, delay_assigned(user_count));
        
        % Uplink delay - first append zeros and then the signal
        Tx_signal_delayed = [delay_zeroes tx_signal_FO(:,1:end-delay_assigned(user_count))];
        
        % Collecting all users signals
        Tx_signal_allUsers(user_count,:) = Tx_signal_delayed;
        
    end
    
    % Summing the signal from all users
    Tx_allUsers_withoutNoise = sum(Tx_signal_allUsers,1);
    
    % SNR loop
    for SNR = SNR_start_dB:SNR_step_dB:SNR_end_dB
        
        SNR_lin = 10^(SNR/10);      % linear value for SNR
        
        % Noise standard deviation
        if(para.pre_domain == "FD")
            noise_stdev = sqrt(signal_power/(2 * SNR_lin * para.N_point))*sqrt(para.num_preamble_sc/para.N_point)*sqrt(para.N_point/(para.N_point+para.cp_length));
        else
            noise_stdev = sqrt(signal_power/(2 * SNR_lin));
        end
        
        % Generating noise for all users and adding together
        noise_dim = size(Tx_allUsers_withoutNoise);
        if(add_noise == 1)
            for kk = 1: para.users
                awgn_noise(kk,:) = (noise_mean + noise_stdev.*randn(noise_dim(1),noise_dim(2)))+1i*(noise_mean + noise_stdev.*randn(noise_dim(1),noise_dim(2)));
            end
        else
            awgn_noise = 0;
        end
        % Summing noises of all users
        awgn_noise_allUsers = sum(awgn_noise,1);
        
        % Finding empirical SNR value
        sig_pow = norm(Tx_allUsers_withoutNoise)^2;
        noise_pow = norm(awgn_noise_allUsers)^2;
        emp_snr(count,SNR_count) = 10*log10(sig_pow/noise_pow);

        
        % Arranging premable(and subband) selected for reference
        for jj = 1: para.num_choices
            [slot_pream{jj},pream_idx{jj}] = sort(pream_select(find(subband_select == jj)));
            [delay_assgn_sort{jj}] = (delay_assigned(find(subband_select == jj))); % This will sort the delay values in ascending order
            delay_assgn_1{jj} = delay_assgn_sort{jj}(pream_idx{jj}); % proper mapping of user count, delay, sub band
            delay_assgn_chanDelay{jj} = delay_assgn_1{jj} + chan_delay;  % Propely arranged delay of each user in the respective sub band
            if( isempty(slot_pream{jj}) )
                slot_pream{jj} = 0;
                delay_assgn_chanDelay{jj} = 0;
                delay_assgn_1{jj} = 0;
            end
        end
        
        
        %% RECEIVER
        
        % Sum of all users' signal and noise at the receiver
        Rx_signal = Tx_allUsers_withoutNoise + awgn_noise_allUsers;
        
        % Considering first n samples of convolution
        Rx_signal = Rx_signal(1,1:numel(tx_signal));
        
        
        % Reshaping to symbol-wise
        if (para.pre_domain == "FD")
            Rx_frame = reshape(Rx_signal,time_grid_size(1),time_grid_size(2));
            % CP removal and FFT
            [frame_dec] = fft_sym_wise(Rx_frame, para); 
           
            [rx_pre_grid] = rx_grid_form(frame_dec,para);
        else
            % Assuming preamble occupies the first symbols
            rx_pre_grid = Rx_signal(1:length(para.para.per_user_pream_seq)).';
        end
        
        
        % Detection
        [vals_1,delay] = preamble_detection(rx_pre_grid, preambles_set, para);
        if(para.pre_domain == "TD")
            vals_1 = vals_1';
        end
        
        
        % Check for all threshold values
        for thresh_ind = 1:1:length(para.threshold)
            
            % det_values gives the detected preambles
            % Considering all values which are more than the threshold as
            % detected preambles
            for ii = 1:1:para.num_choices
                det_values{ii} = (find(vals_1(:,ii) > para.threshold(thresh_ind)))';
                if( isempty(det_values{ii}) )
                    det_values{ii} = 0;
                end
            end
            
            % Checking for false positives
            % False positive is when atleast one detected value isn't a chosen value
            for ii = 1:1:para.num_choices
                temp1 = det_values{ii};
                temp2 = slot_pream{ii};
                fa_flag(ii) = 0;
                if (det_values{ii}~=0)
                    if (sum(ismember(temp1,temp2))~=length(temp1))
                        %false_pos(thresh_ind,SNR_count,ii) = false_pos(thresh_ind,SNR_count,ii) + 1;
                        false_pos(SNR_count,ii,thresh_ind) = false_pos(SNR_count,ii,thresh_ind) + 1;
                        fa_flag(ii) = 1;
                        
                        % Removing the false preamble values from the
                        % detected set
                        det_values{ii} = det_values{ii}(find(ismember(det_values{ii},slot_pream{ii})));
                    end
                end
            end
            
            % Checking for collisions when preambles are selected
            % Collision is when 2 or more users choose the same preamble
            for ii = 1:1:para.num_choices
                if (length(unique(slot_pream{ii}))<length(slot_pream{ii}))
                    coll(SNR_count,ii,thresh_ind) = coll(SNR_count,ii,thresh_ind) + 1;
                    %coll(thresh_ind,SNR_count,ii) = coll(thresh_ind,SNR_count,ii) + 1;
                end
            end
            
            % Checking if there is collision and then decide all user succ or not
            % Checking if the delays of detected preambles are within
            % tolerance
            succ_SB = zeros(1,para.num_choices);
            succ_delayTol_SB = zeros(1,para.num_choices);
            for ii = 1:1:para.num_choices
                if (det_values{ii}~=0)
                    if(length(unique(slot_pream{ii}))==length(slot_pream{ii}))
                        % Checking if chosen and det preambles set match exactly
                        % If they match, all users are successful
                        if(isequal(slot_pream{ii},det_values{ii}))
                            
                            if(para.pre_domain == "TD")
                                delay_values = (delay(det_values{ii}));
                                % Delay values of all users
                                delay_allUsers = delay_values - para.pre_length;
                            else
                                delay_values = (delay(det_values{ii},ii))';   % Detected delay values
                                delay_allUsers = delay_values;
                            end
                            
                            % delay_diff is the difference between delay values detected and delay values assigned to all users
                            delay_diff = delay_allUsers - delay_assgn_chanDelay{ii};
                            
                            % Taking count of all users successful, for taking probability of succ later
                            succ_SB(ii) = 1;
                            
                            % Checking with delay tolerance
                            if(abs(delay_diff) <= delay_tolerance)
                                succ_delayTol_SB(ii) = 1;
                            else
                                disp('delay is more');
                                %[vals_1,delay] = preamble_detection(rx_pre_grid, preambles_set, para);
                            end
                            
                            %Storing the delay value
                            store_delay_diff{ii} =  delay_diff;
                        end
                    end
                elseif (slot_pream{ii}==0)
                    % If no preamble is chosen in a sub band and nothing is
                    % detected, then putting the value as 1 will make sure
                    % it is right. The logic is if proper detection, then 1,
                    % else 0. So if after all sub bands results are there,
                    % total should add up to number of subbands, then all
                    % users are successful
                    succ_SB(ii) = 1;    
                    succ_delayTol_SB(ii) = 1;
                end
                
            end
            
            % Storing the delay difference
            delay_diff_list{thresh_ind}{count,SNR_count} = cell2mat(store_delay_diff);
            
            % If in all sub bands, users are successful, then this count is
            % taken as all users are successful
            if(sum(fa_flag)==0)
                if(sum(succ_delayTol_SB)==para.num_choices)
                    count_succ(SNR_count,thresh_ind) = count_succ(SNR_count,thresh_ind) + 1;
                end
            end
             
        end
        
        
        %Collecting false positives
        false_pos_table = false_pos;
        
        %Collecting the count of all users being successful
        count_succ_table = count_succ;
        
        %Collecting collisions
        coll_table = coll;
        
        SNR_count = SNR_count + 1;
        
    end
    
    % Saving the values in a mat file
    save('TD_lteFC.mat','SNR_range','false_pos_table','count_succ_table','coll_table','delay_diff_list','emp_snr');
toc
end

% Empirical SNR value
emp_SNR_list = mean(emp_snr);

% Probability of all users being successful
prob_count_succ_table = count_succ./no_iter;

% Probability of false positives
prob_false_pos_table = false_pos./no_iter;

% Collisions probability
prob_coll_table = coll./no_iter;

% Prob of false pos and collision - summed up along all subbands - for
% plotting purpose
for ii = 1:1:length(para.threshold)
    prob_FP(:,ii) = sum(prob_false_pos_table(:,:,ii),2);
    prob_coll(:,ii) = sum(prob_coll_table(:,:,ii),2); 
end


 

%% Plots

%CDF plot for the delay difference

for ii = 1:1:length(para.threshold)
    delay_diff_eachThresh = delay_diff_list{ii};
    delay_diff_allValues = [];
    for jj = 1:1:no_iter
        delay_diff_allValues = [cell2mat(delay_diff_eachThresh(jj,:)) delay_diff_allValues];
    end
    cdf_thresh = cdfplot(delay_diff_allValues);
    hold on;
    ylim([-0.5 1.5]);
  
    legend('Threshold=3','Threshold = 3.1','Threshold=3.2','Threshold = 3.3','Threshold = 3.4','Threshold = 3.5','Threshold = 3.6','Threshold=3.8','Threshold = 4','Threshold=4.2','Threshold=4.4','Threshold=4.5','Threshold=4.6','Threshold = 4.8','Threshold =5','Threshold=5.2 ');
    %legend('Threshold = 80');
    set(cdf_thresh,'MarkerSize', 10, 'LineWidth', 2);
end



% Plotting probabilities - false pos and all user succ

x = emp_SNR_list;
%if(add_noise==0)
%    x=[SNR_start_dB: SNR_step_dB:SNR_end_dB];
%end
%Plot of SNR v/s probability of false positives
figure(2);
semilogy(x, prob_FP, '+-', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('SNR per user in dB');
ylabel('Prob of false positives');
title('To find threshold-FD');
%legend('Threshold = 56','Threshold=57', 'Threshold=58', 'Threshold=59', 'Threshold=60');
legend('Threshold = 5','Threshold=6', 'Threshold=7', 'Threshold=8', 'Threshold=9', 'Threshold=10');
%legend('Threshold=3','Threshold = 3.1','Threshold=3.2','Threshold = 3.3','Threshold = 3.4','Threshold = 3.5','Threshold = 3.6','Threshold=3.8','Threshold = 4','Threshold=4.2','Threshold=4.4','Threshold=4.5','Threshold=4.6','Threshold = 4.8','Threshold =5','Threshold=5.2 ');

figure(3);
semilogy(x, prob_count_succ_table, '+-', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
xlabel('SNR per user in dB');
ylabel('Prob of all users detected');
title('Prob of all users successful-within delay tolerance-FD');
legend('Threshold = 5','Threshold=6', 'Threshold=7', 'Threshold=8', 'Threshold=9', 'Threshold=10');
%legend('Threshold = 56','Threshold=57', 'Threshold=58', 'Threshold=59', 'Threshold=60');
%legend('Threshold=3','Threshold = 3.1','Threshold=3.2','Threshold = 3.3','Threshold = 3.4','Threshold = 3.5','Threshold = 3.6','Threshold=3.8','Threshold = 4','Threshold=4.2','Threshold=4.4','Threshold=4.5','Threshold=4.6','Threshold = 4.8','Threshold =5','Threshold=5.2 ');

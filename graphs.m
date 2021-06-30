data1=load('sb6_long_all_thresholds_both');
data2=load('sb6_multi_no_spread_all_thresholds_both');

%for ii=1:1:6
figure(2);
semilogy(data1.x, data1.prob_FP(:,6), '+-', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
semilogy(data2.x, data2.prob_FP(:,7), '+-', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
%end
grid on;
xlabel('SNR per user in dB');
ylabel('Prob of false positives');
%title('Prob of all users successful-within delay tolerance and only noise-FD');
%title('Prob of all users successful-within delay tolerance and noise and channel-FD');
%title('1 subband choice, 1000 iterations with only noise and threshold=3.8');
title('1 subband choices, 1000 iterations with both noise and channel and threshold=3.8');
%legend('Threshold = 80');
legend('long','multiple short without spread sequence','multiple short spread sequence');
%legend('Threshold=3','Threshold = 3.1','Threshold=3.2','Threshold = 3.3','Threshold = 3.4','Threshold = 3.5','Threshold = 3.6','Threshold=3.8');
%legend('Threshold=3','Threshold = 3.1','Threshold=3.2','Threshold = 3.3','Threshold = 3.4','Threshold = 3.5','Threshold = 3.6','Threshold=3.8','Threshold = 4','Threshold=4.2','Threshold=4.4','Threshold=4.5','Threshold=4.6','Threshold = 4.8','Threshold =5','Threshold=5.2 ');
figure(3);
semilogy(data1.x, data1.prob_count_succ_table(:,6), '+-', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
semilogy(data2.x, data2.prob_count_succ_table(:,7), '+-', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
%end
grid on;
xlabel('SNR per user in dB');
ylabel('Prob of all users detected');
%title('Prob of all users successful-within delay tolerance and only noise-FD');
%title('Prob of all users successful-within delay tolerance and noise and channel-FD');
%title('1 subband choice, 1000 iterations with only noise and threshold=3.8');
%title('1 subband choice, 1000 iterations with both noise and channel and threshold=3.8');
%legend('Threshold = 80');
title('comparison between long and multiple short for subband choices=6');
legend('long preamble threshold=3.5','multiple-short preamble threshold=3.6');
%legend('long','multiple short without spread','multiple short spread_sequence');
%legend('Threshold=3','Threshold = 3.1','Threshold=3.2','Threshold = 3.3','Threshold = 3.4','Threshold = 3.5','Threshold = 3.6','Threshold=3.8');
%legend('Threshold=3','Threshold = 3.1','Threshold=3.2','Threshold = 3.3','Threshold = 3.4','Threshold = 3.5','Threshold = 3.6','Threshold=3.8','Threshold = 4','Threshold=4.2','Threshold=4.4','Threshold=4.5','Threshold=4.6','Threshold = 4.8','Threshold =5','Threshold=5.2 ');


%for ii=1:1:6
figure(4);
for ii=1:1:10
    if ii~=7
    semilogy(data2.x, data2.prob_FP(:,ii), '+-', 'MarkerSize', 10, 'LineWidth', 2);
    end
%semilogy(data2.x, data1.prob_count_succ_table(:,ii), '+-', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
%semilogy(data2.x, data2.prob_count_succ_table(:,8), '+-', 'MarkerSize', 10, 'LineWidth', 2);
%hold on;
end
grid on;
xlabel('SNR per user in dB');
ylabel('Prob of all users detected');
%title('Prob of all users successful-within delay tolerance and only noise-FD');
%title('Prob of all users successful-within delay tolerance and noise and channel-FD');
%title('1 subband choice, 1000 iterations with only noise and threshold=3.8');
title('To find threshold-FD');
%title('1 subband choice, 1000 iterations with both noise and channel and threshold=3.8');
%legend('Threshold = 80');
%legend('long','multiple short without spread','multiple short spread_sequence');
%legend('Threshold = 3.6','Threshold=3.8','Threshold=4');
legend('Threshold=3','Threshold = 3.1','Threshold=3.2','Threshold = 3.3','Threshold = 3.4','Threshold = 3.5','Threshold = 3.6','Threshold=3.8','Threshold = 4','Threshold=4.2','Threshold=4.4','Threshold=4.5','Threshold=4.6','Threshold = 4.8','Threshold =5','Threshold=5.2 ');

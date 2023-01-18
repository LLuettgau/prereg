close all
restoredefaultpath

%% --------------------------------------------------------------------------
% Load pre-computed effects and define samples
save_plots = 0;

PRE = load([is.OutputDir 'resting_results/', 'StimAll_rst0']);
POST = load([is.OutputDir 'resting_results/', 'StimAll_rst1']);

% Final sample
conId = 1;
included = [1 3 4 5 6 7 8 9 10 11]; %which subjects to include, exclude subj2

samples = {};
samples{1,1} = conId; samples{1,2} = 'b';

% FWD - BWD [SEQUENCENESS]
PRE_seq =  PRE.sfAll(2:end, :, included,:,:,:) -  PRE.sbAll(2:end, :, included,:,:,:); % 1st lag is 0ms (NaN)
POST_seq = POST.sfAll(2:end, :, included,:,:,:) - POST.sbAll(2:end, :, included,:,:,:);

PRE_POST_seq_diff = POST_seq - PRE_seq;

%prepare for plotting p-values
value_grid = NaN(size(PRE_POST_seq_diff,5),size(PRE_POST_seq_diff,6),3,2);
%% --------------------------------------------------------------------------
% Sequenceness x lag x training time x L1 regularization params
for iTT = 1:size(PRE_POST_seq_diff,5)    
    for l1reg = 1:size(PRE_POST_seq_diff,6)
        if l1reg == 1
            counter = 1;
            h = figure;
        elseif l1reg == (size(PRE_POST_seq_diff,6) / 2) + 1 
            counter = 1;
            h = figure;
        end
        
        set(gcf, 'Units', 'point', 'Position', [50 50 700, 1000])
        for isess = 2
            if isess == 1
                b = PRE_seq(:,:,:,:,iTT,l1reg); % [lag, perm, subj, irrelevant, training time, L1 param]
                tt = 'PRE-learning rest';
            elseif isess == 2
                b =  POST_seq(:,:,:,:,iTT,l1reg);
                tt = 'MID-learning rest';
            elseif isess == 3
                b =  PRE_POST_seq_diff(:,:,:,:,iTT,l1reg);
                tt = 'PRE-MID-difference';
            end
                        
            mn = mean(b(:, 1,:), 3);
            mean_null = mean(mean(b(:, 2:end, :),3),2);
            
            sem = std(squeeze(b(:, 1, :)),0,2)./sqrt(length(squeeze(b(:, 1, :))));
            
            shadedErrorBar(10:10:600, mn, sem, samples{1, 2}, 1);
            
            % for each permutation, calculate mean absolute sequenceness over participants,
            % then take the maximal value over all perms and lags
            hold on;
            
%             thresh_lag = prctile(abs(mean(b(:, 2:end, :),3)),95,2);
            thresh_lag = prctile(mean(b(:, 2:end, :),3),95,2);

            thresh = max(thresh_lag,[],1);
            
            plot([0 600], [thresh, thresh], [samples{1,2} ':'], 'LineWidth', 2)
            
            plot(10:10:600, thresh_lag, [samples{1,2} '-'], 'LineWidth', 2)
           
            
            %calculate and save p-value of sequenceness effect at 40-60ms in mid rest period
            if isess == 2
                %count the number of values (statistics) that are greater than or equal to the observed value,
                %and divide by the number of values. pval = sum(s >= s0)/N;
                null_means = mean(b([4 5 6], 2:end, :), 3);
                null_dist = max(null_means, [], 1);
                figure, hist(null_dist)
                realval_40ms = mean(b(4, 1, :));
                realval_50ms = mean(b(5, 1, :));
                realval_60ms = mean(b(6, 1, :));

                PVALS(1) = mean(realval_40ms <= null_dist);
                PVALS(2) = mean(realval_50ms <= null_dist);
                PVALS(3) = mean(realval_60ms <= null_dist);
                FINAL_PVAL = min(PVALS);
                
                hold on, plot(realval_40ms * [1 1], [0 250])
                % print empirical p-values of sequenceness effect (one-tailed)
                disp(["FINAL P VALUE = " num2str(FINAL_PVAL)])
             end
            
            if isess == 1
               ylabel('sequenceness', 'FontSize', 16)
            end

            
            if l1reg == 1 || l1reg == (size(PRE_POST_seq_diff,6) / 2) + 1  
                title([tt ', training time: ' num2str(is.tss(iTT) * 10) ' ms'], 'FontSize', 16)
            end
           
            if l1reg == size(PRE_POST_seq_diff,6)
                xlabel('lag (ms)', 'FontSize', 16)
            end
            screen_size = get(0, 'ScreenSize');
            
        end
        
    end
end




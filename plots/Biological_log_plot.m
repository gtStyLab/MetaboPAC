% Plot results when using response factors sampled from a log uniform
% distribution on noiseless biological systems.
% Figure S10

clear;clc;

analysis_count = 1;
percKnownKinetics = 100;
for run = 1:20
    fileName = sprintf('results/yeast_MetaboPAC_log_percKnownKinetics-%03d_rep-%03d.mat',percKnownKinetics,run);

    if isfile(fileName)
        load(fileName,'RF_opt','trueRF','numMetabs','maxRandVal','minRandVal','knownMet');
        
        std_F(analysis_count,:) = std(RF_opt(:,1:numMetabs));
        med_F(analysis_count,:) = median(RF_opt(:,1:numMetabs));
        trueRF_list(analysis_count,:) = trueRF;
        %rand_F(analysis_count,:) = (maxRandVal-minRandVal).*rand(1,numMetabs) + minRandVal;
        trueRF_lower_mag = 0;
        trueRF_upper_mag = 3;
        rand_F(analysis_count,:) = 10.^(trueRF_lower_mag + (trueRF_upper_mag-trueRF_lower_mag)*rand(1,numMetabs));
        F500(analysis_count,:) = 500*ones(1,numMetabs);
        
        abs_log2_error(analysis_count,:) = abs(log2(med_F(analysis_count,:)) - log2(trueRF));
        abs_log2_error_random(analysis_count,:) = abs(log2(rand_F(analysis_count,:)) - log2(trueRF));
        abs_log2_error_F500(analysis_count,:) = abs(log2(F500(analysis_count,:)) - log2(trueRF));
        
        analysis_count = analysis_count + 1;
    end
end

% Log2 error plot
f = figure(1);
f.Position = [500 500 1400 600];
annotation('line', [0.49 0.49], [0.05, 0.95])
a = annotation('textbox', [0.13 0.99 0 0],'String',{'A.'},'EdgeColor','none');
a.FontSize = 16;
a.FontWeight = 'bold';
b = annotation('textbox', [0.57 0.99 0 0],'String',{'B.'},'EdgeColor','none');
b.FontSize = 16;
b.FontWeight = 'bold';
subplot(1,2,1);
y = 100*[mean(sum(abs_log2_error < log2(1.1),2)/numMetabs) mean(sum(abs_log2_error_random < log2(1.1),2)/numMetabs) mean(sum(abs_log2_error_F500 < log2(1.1),2)/numMetabs);...
    mean(sum(abs_log2_error < log2(1.3),2)/numMetabs) mean(sum(abs_log2_error_random < log2(1.3),2)/numMetabs) mean(sum(abs_log2_error_F500 < log2(1.3),2)/numMetabs);...
    mean(sum(abs_log2_error < log2(1.5),2)/numMetabs) mean(sum(abs_log2_error_random < log2(1.5),2)/numMetabs) mean(sum(abs_log2_error_F500 < log2(1.5),2)/numMetabs)];
bar(y);
set(gca, 'XTickLabel',{'< log_2(1.1)','< log_2(1.3)','< log_2(1.5)'});
hold on

% Error bars
err = 100*[std(sum(abs_log2_error < log2(1.1),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_random < log2(1.1),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_F500 < log2(1.1),2)/numMetabs)/sqrt(analysis_count-1);...
    std(sum(abs_log2_error < log2(1.3),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_random < log2(1.3),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_F500 < log2(1.3),2)/numMetabs)/sqrt(analysis_count-1);...
    std(sum(abs_log2_error < log2(1.5),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_random < log2(1.5),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_F500 < log2(1.5),2)/numMetabs)/sqrt(analysis_count-1)];
ngroups = 3;
nbars = 3;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), err(:,i), 'k.');
end
ylim([0 105])


legend('MetaboPAC','Random','500','Location','Northwest');

title({'S. cerevisiae',sprintf('%d%% of kinetic equations known',percKnownKinetics)})
%ylabel({'Percent of predicted conversion';'factors within \chi% error'});
ylabel({'% of predicted response','factors within error range'});
set(gca,'FontSize',15)


clear;clc;

analysis_count = 1;
percKnownKinetics = 100;
for run = 1:20
    fileName = sprintf('results/ecoli_MetaboPAC_log_percKnownKinetics-%03d_rep-%03d.mat',percKnownKinetics,run);

    if isfile(fileName)
        load(fileName,'RF_opt','trueRF','numMetabs','maxRandVal','minRandVal','knownMet');
        
        std_F(analysis_count,:) = std(RF_opt(:,1:numMetabs));
        med_F(analysis_count,:) = median(RF_opt(:,1:numMetabs));
        trueRF_list(analysis_count,:) = trueRF;
        %rand_F(analysis_count,:) = (maxRandVal-minRandVal).*rand(1,numMetabs) + minRandVal;
        trueRF_lower_mag = 0;
        trueRF_upper_mag = 3;
        rand_F(analysis_count,:) = 10.^(trueRF_lower_mag + (trueRF_upper_mag-trueRF_lower_mag)*rand(1,numMetabs));
        F500(analysis_count,:) = 500*ones(1,numMetabs);
        
        abs_log2_error(analysis_count,:) = abs(log2(med_F(analysis_count,:)) - log2(trueRF));
        abs_log2_error_random(analysis_count,:) = abs(log2(rand_F(analysis_count,:)) - log2(trueRF));
        abs_log2_error_F500(analysis_count,:) = abs(log2(F500(analysis_count,:)) - log2(trueRF));
        
        analysis_count = analysis_count + 1;
    end
end

% Log2 error plot
figure(1);
subplot(1,2,2);
y = 100*[mean(sum(abs_log2_error < log2(1.1),2)/numMetabs) mean(sum(abs_log2_error_random < log2(1.1),2)/numMetabs) mean(sum(abs_log2_error_F500 < log2(1.1),2)/numMetabs);...
    mean(sum(abs_log2_error < log2(1.3),2)/numMetabs) mean(sum(abs_log2_error_random < log2(1.3),2)/numMetabs) mean(sum(abs_log2_error_F500 < log2(1.3),2)/numMetabs);...
    mean(sum(abs_log2_error < log2(1.5),2)/numMetabs) mean(sum(abs_log2_error_random < log2(1.5),2)/numMetabs) mean(sum(abs_log2_error_F500 < log2(1.5),2)/numMetabs)];
bar(y);
set(gca, 'XTickLabel',{'< log_2(1.1)','< log_2(1.3)','< log_2(1.5)'});
hold on

% Error bars
err = 100*[std(sum(abs_log2_error < log2(1.1),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_random < log2(1.1),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_F500 < log2(1.1),2)/numMetabs)/sqrt(analysis_count-1);...
    std(sum(abs_log2_error < log2(1.3),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_random < log2(1.3),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_F500 < log2(1.3),2)/numMetabs)/sqrt(analysis_count-1);...
    std(sum(abs_log2_error < log2(1.5),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_random < log2(1.5),2)/numMetabs)/sqrt(analysis_count-1) std(sum(abs_log2_error_F500 < log2(1.5),2)/numMetabs)/sqrt(analysis_count-1)];
ngroups = 3;
nbars = 3;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), err(:,i), 'k.');
end
ylim([0 105])


title({'E. coli',sprintf('%d%% of kinetic equations known',percKnownKinetics)})
%ylabel({'Percent of predicted conversion';'factors within \chi% error'});
ylabel({'% of predicted response','factors within error range'});
set(gca,'FontSize',15)



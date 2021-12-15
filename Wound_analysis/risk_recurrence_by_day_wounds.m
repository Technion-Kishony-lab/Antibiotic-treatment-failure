function [] = risk_recurrence_by_day_wounds(wound_cases,params)
%% sensitive chance of failure
% function takes wound_cases structure and optional params and makes plot of 
% the 7-day moving avaerage % of early recurrences for all antibtioic treated cases   

wound_cases.hasdiag = true(size(wound_cases.RandomID));
av_window = 7;
days = 4:1:50;
%for all days
for kk = 1:length(days)-1
clear treatfailure
treatfailure = wound_cases.date_diff >= days(kk) & wound_cases.date_diff < days(kk)+av_window;
%for all antibioitcs
for drug = 1:params.number_drugs
 
   all_sensitive1 = ismember(wound_cases.SMP_Res(:,drug), params.sensitive_group) & wound_cases.hasdiag;
   all_resistant1 = ismember(wound_cases.SMP_Res(:,drug), params.resistant_group) & wound_cases.hasdiag;
   all_sensitive_nextres1 =  ismember(wound_cases.next_res(:,drug), params.sensitive_group) ;
   all_resistant_nextres1 =  ismember(wound_cases.next_res(:,drug), params.resistant_group) ;
   all_nextres = ismember(wound_cases.next_res(:,drug),[1 2 3]);
   
   all_sensitive_purchased1(drug ,kk) =  nnz((all_sensitive1 | all_resistant1) & wound_cases.PCR_sameday(:,drug));
   sensitive_purchased_fails1(drug ,kk) = nnz(all_sensitive1 & wound_cases.PCR_sameday(:,drug) & treatfailure);
   sensitive_purchased_nextresfails1(drug ,kk) = nnz(all_sensitive1 & wound_cases.PCR_sameday(:,drug) & treatfailure & all_nextres);
   sensitive_purchased_fails_SS(drug ,kk) = nnz(all_sensitive1 & wound_cases.PCR_sameday(:,drug) & treatfailure & all_sensitive_nextres1);
   sensitive_purchased_fails_SR(drug ,kk) = nnz(all_sensitive1 & wound_cases.PCR_sameday(:,drug) & treatfailure & all_resistant_nextres1);
     
   all_resistant_purchased(drug ,kk) = nnz((all_sensitive1 | all_resistant1) & wound_cases.PCR_sameday(:,drug));
   resistant_purchased_fails(drug ,kk) = nnz(all_resistant1 & wound_cases.PCR_sameday(:,drug) & treatfailure);
   resistant_purchased_nextresfails(drug ,kk) = nnz(all_resistant1 & wound_cases.PCR_sameday(:,drug) & treatfailure & all_nextres);
   resistant_purchased_fails_RS(drug ,kk) = nnz(all_resistant1 & wound_cases.PCR_sameday(:,drug) & treatfailure & all_sensitive_nextres1);
   resistant_purchased_fails_RR(drug ,kk) = nnz(all_resistant1 & wound_cases.PCR_sameday(:,drug) & treatfailure & all_resistant_nextres1);
   
    
end

% sum the numbers from all the drugs for each day
fail_split1_all(1:2,kk, 1) = [sum(sensitive_purchased_fails_SS(1:params.number_drugs ,kk))/sum(sensitive_purchased_nextresfails1(1:params.number_drugs ,kk)) ; sum(sensitive_purchased_fails_SR(1:params.number_drugs ,kk))/sum(sensitive_purchased_nextresfails1(1:params.number_drugs ,kk))]*(sum(sensitive_purchased_fails1(1:params.number_drugs ,kk))/sum(all_sensitive_purchased1(1:params.number_drugs))*100);
fail_split_res1_all(1:2,kk, 1) = [sum(resistant_purchased_fails_RR(1:params.number_drugs ,kk))/sum(resistant_purchased_nextresfails(1:params.number_drugs ,kk)) ; sum(resistant_purchased_fails_RS(1:params.number_drugs ,kk))/sum(resistant_purchased_nextresfails(1:params.number_drugs ,kk))]*(sum(resistant_purchased_fails(1:params.number_drugs ,kk))/sum(all_resistant_purchased(1:params.number_drugs))*100);

      
end

%% plot the fig
fail_split2(:,:,1) = fail_split1_all(:,:,1);
fail_split_res2(:,:,1) = fail_split_res1_all(:,:,1);

b1 = bar(1:length(days)-1,[fail_split2(:,:,1)',fail_split_res2(:,:,1)']./av_window, 'stacked' , 'FaceColor','flat',  'EdgeColor', 'none', 'BarWidth', 1);

b1(1).CData(:,:) = ones(length(days)-1,1)*params.SS_color;
b1(2).CData(:,:) = ones(length(days)-1,1)*params.SR_color;
b1(3).CData(:,:) = ones(length(days)-1,1)*params.RR_color;
b1(4).CData(:,:) = ones(length(days)-1,1)*params.RS_color;

ylabel('% of patients');
set(gca,'XTick',1:8:length(days));
days_string = num2str((days(1:8:length(days)))');
set(gca,'xticklabel',days_string);
xlabel('# of days between first and second sample');
xlim([0.5 length(days)-0.5])
hold on
ylim([0 0.3]);
yl = ylim;
plot([28-days(1)+1 28-days(1)+1],yl, 'k--');
end
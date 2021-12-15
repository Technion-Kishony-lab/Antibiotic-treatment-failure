function []= mismatched_vs_matched_recurrence_rate(UTI_cases,params)
% function takes UTI_case struc and optional params and makes bar chart of 
% the % of early recurrences for sensitivity matched and mismatched antibiotic 
% treated cases.

% pre-allocate
total_prch = zeros(1,params.number_drugs);
all_sensitive_purchased = zeros(1,params.number_drugs);
sensitive_purchased_fails = zeros(1,params.number_drugs);
sensitive_purchased_nextresfails = zeros(1,params.number_drugs);
sensitive_purchased_fails_SS = zeros(1,params.number_drugs);
sensitive_purchased_fails_SR = zeros(1,params.number_drugs);
all_resistant_purchased = zeros(1,params.number_drugs);
resistant_purchased_fails = zeros(1,params.number_drugs);
resistant_purchased_nextresfails = zeros(1,params.number_drugs);
resistant_purchased_fails_RS = zeros(1,params.number_drugs);
resistant_purchased_fails_RR = zeros(1,params.number_drugs);

% loop over all antibiotics
for drug = 1:params.number_drugs  
   total_prch(drug) = nnz(UTI_cases.PCR_sameday(:,drug)  & UTI_cases.hasdiag );
   all_sensitive = (ismember(UTI_cases.SMP_Res(:,drug), params.sensitive_group)) & UTI_cases.hasdiag  ;
   all_reistant = (ismember(UTI_cases.SMP_Res(:,drug), params.resistant_group))  & UTI_cases.hasdiag ;
      
   all_sensitive_nextres =  ismember(UTI_cases.next_res(:,drug),params.sensitive_group);
   all_resistant_nextres =  ismember(UTI_cases.next_res(:,drug),params.resistant_group);
   all_nextres = ismember(UTI_cases.next_res(:,drug),[1 2 3]);
   
   all_sensitive_purchased(drug) =  nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug));
   sensitive_purchased_fails(drug) = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure);
   sensitive_purchased_nextresfails(drug) = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_nextres);
   
   sensitive_purchased_fails_SS(drug) = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_sensitive_nextres);
   sensitive_purchased_fails_SR(drug) = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_resistant_nextres);
   
   all_resistant_purchased(drug) =  nnz(all_reistant & UTI_cases.PCR_sameday(:,drug));
   resistant_purchased_fails(drug) = nnz(all_reistant & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure);
   resistant_purchased_nextresfails(drug) = nnz(all_reistant & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_nextres);
   
   resistant_purchased_fails_RS(drug) = nnz(all_reistant & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_sensitive_nextres);
   resistant_purchased_fails_RR(drug) = nnz(all_reistant & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_resistant_nextres);
end

% sum the numbers from all antibiotics 
all_ratio = sum(sensitive_purchased_fails)/sum(all_sensitive_purchased)*100;
all_fail_split = [sum(sensitive_purchased_fails_SS)/sum(sensitive_purchased_nextresfails) ; sum(sensitive_purchased_fails_SR)/sum(sensitive_purchased_nextresfails)]*all_ratio;
all_fail_split_error = sqrt((all_ratio/100).*(1-(all_ratio/100))./(sum(all_sensitive_purchased)))*100;
 
all_resistant_ratio = sum(resistant_purchased_fails)/sum(all_resistant_purchased)*100;
all_resistant_fail_split = [sum(resistant_purchased_fails_RS)/sum(resistant_purchased_nextresfails) ; sum(resistant_purchased_fails_RR)/sum(resistant_purchased_nextresfails)]*all_resistant_ratio;
all_resistant_fail_split_error = sqrt((all_resistant_ratio/100).*(1-(all_resistant_ratio/100))./(sum(all_resistant_purchased)))*100;

%% print bar chart 
figure
HW = [3 2];
set(gcf,'color','w')
clf; set(gcf,'name','Fig. 1B');
set(gcf,'units','centimeters','Position',[2 2 7 10]);

width = 0.5;

ii = 1;
    hold on
xpos = 2;%gap+0.32;
hd1 = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],...
    [0  0 all_fail_split(2,ii) all_fail_split(2,ii)],...
    'w','EdgeColor','k', 'FaceColor', params.SR_color); 
hd2 = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],...
    [all_fail_split(2,ii) all_fail_split(2,ii) ...
    (all_fail_split(1,ii)+all_fail_split(2,ii)) ...
    (all_fail_split(1,ii)+all_fail_split(2,ii))],...
    'w','EdgeColor','k', 'FaceColor', params.SS_color); 
errorbar(xpos,(all_fail_split(1,ii)+all_fail_split(2,ii)),all_fail_split_error(ii), 'k', 'LineStyle','none')

hold on
xpos = 1;%gap-0.32;
hd1 = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],...
    [0  0 all_resistant_fail_split(2,ii) all_resistant_fail_split(2,ii)],...
    'w','EdgeColor','k', 'FaceColor', params.RR_color); 
hd2 = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],...
    [all_resistant_fail_split(2,ii) all_resistant_fail_split(2,ii) ...
    (all_resistant_fail_split(1,ii)+all_resistant_fail_split(2,ii)) ...
    (all_resistant_fail_split(1,ii)+all_resistant_fail_split(2,ii))],...
    'w','EdgeColor','k', 'FaceColor', params.RS_color); 
errorbar(xpos,(all_resistant_fail_split(1,ii)+all_resistant_fail_split(2,ii)),all_resistant_fail_split_error(ii), 'k', 'LineStyle','none')

ylabel({'percentage of UTIs' ; 'resulting in early recurrence'})
set(gca,'XTick',[1 2]);
set(gca,'xticklabel',{'mismatched','mathced'});
xtickangle(45);

ylim([0 25])
end
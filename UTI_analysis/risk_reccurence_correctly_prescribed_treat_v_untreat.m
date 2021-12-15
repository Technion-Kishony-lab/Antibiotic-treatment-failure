function []= risk_reccurence_correctly_prescribed_treat_v_untreat(UTI_cases,params)
% function takes UTI_case structure and optional params and makes bar charts of 
% the % of early recurrences sensitive to specific antiboitic, for untreated cases and 
% cases treated the the specific antibtioic  

%loop all antibiotics
for drug = 1:params.number_drugs  
 
   total_prch(drug) = nnz(UTI_cases.PCR_sameday(:,drug) & UTI_cases.hasdiag);
   all_sensitive = (ismember(UTI_cases.SMP_Res(:,drug), params.sensitive_group)) & UTI_cases.hasdiag ;   
   all_sensitive_nextres =  ismember(UTI_cases.next_res(:,drug),params.sensitive_group);
   all_resistant_nextres =  ismember(UTI_cases.next_res(:,drug),params.resistant_group);
   all_nextres = ismember(UTI_cases.next_res(:,drug),[1 2 3]);
   
   all_sensitive_purchased(drug) =  nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug));
   sensitive_purchased_fails(drug) = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure);
   sensitive_purchased_nextresfails(drug) = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_nextres);
   
   sensitive_purchased_fails_SS(drug) = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_sensitive_nextres);
   sensitive_purchased_fails_SR(drug) = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_resistant_nextres);
   
   %first calulate rate of failure 
   ratio(drug) = sensitive_purchased_fails(drug)/all_sensitive_purchased(drug)*100;
   number_prch(drug) = all_sensitive_purchased(drug);
   ratio_error(drug) = sqrt((ratio(drug)/100).*(1-(ratio(drug)/100))./(all_sensitive_purchased(drug)))*100;
   % slit the overall rate of faliure into 2 modes based on ratio of gained
   % res to remained sen
   fail_split(1:2,drug) = [sensitive_purchased_fails_SS(drug)/sensitive_purchased_nextresfails(drug) ; sensitive_purchased_fails_SR(drug)/sensitive_purchased_nextresfails(drug)]*ratio(drug);
   
   %same for untreated cases  
   all_sensitive_noprch =  nnz(all_sensitive & UTI_cases.PCR_sameday(:,10));
   sensitive_noprch_fails = nnz(all_sensitive & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure );
   sensitive_noprch_nextresfails = nnz(all_sensitive & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure & all_nextres);
   sensitive_noprch_fails_SS = nnz(all_sensitive & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure & all_sensitive_nextres);
   sensitive_noprch_fails_SR = nnz(all_sensitive & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure & all_resistant_nextres);
   
   ratio_noprch(drug) = sensitive_noprch_fails/all_sensitive_noprch*100;
   number_noprch(drug) = all_sensitive_noprch;
   ratio_noprch_error(drug) = sqrt((ratio_noprch(drug)/100).*(1-(ratio_noprch(drug)/100))./(all_sensitive_noprch))*100;

   fail_split_noprch(1:2,drug) = [sensitive_noprch_fails_SS/sensitive_noprch_nextresfails ; sensitive_noprch_fails_SR/sensitive_noprch_nextresfails]*ratio_noprch(drug);
   
end

%% plot figure

width = 0.5;
gap = 1.8;
plot([0 7*gap+3], [0 0],'k')

for jj = 1:params.number_drugs
    ii = params.new_order(jj);
    hold on
xpos = jj*gap-0.32;
hd3(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 fail_split_noprch(2,ii) fail_split_noprch(2,ii)],'w','EdgeColor','k', 'FaceColor', params.SR_color); % using MATLAB "patch" to establish the border
hd4(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[fail_split_noprch(2,ii) fail_split_noprch(2,ii) (fail_split_noprch(1,ii)+fail_split_noprch(2,ii)) (fail_split_noprch(1,ii)+fail_split_noprch(2,ii))],'w','EdgeColor','k', 'FaceColor', params.SS_color); % using MATLAB "patch" to establish the border
errorbar(xpos,(fail_split_noprch(1,ii)+fail_split_noprch(2,ii)),ratio_noprch_error(ii), 'k', 'LineStyle','none')
end

for jj = 1:params.number_drugs
    ii = params.new_order(jj);
hold on
xpos = jj*gap+0.32;
hd1(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 fail_split(2,ii) fail_split(2,ii)],'w','EdgeColor','k', 'FaceColor', params.SR_color); % using MATLAB "patch" to establish the border
hd2(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[fail_split(2,ii)  fail_split(2,ii) (fail_split(1,ii)+fail_split(2,ii)) (fail_split(1,ii)+fail_split(2,ii))],'w','EdgeColor','k', 'FaceColor', params.SS_color); % using MATLAB "patch" to establish the border
hatch(hd3(jj),[45,2.5,0.8],'k')
hatch(hd4(jj),[45,2.5,0.8],'w')
errorbar(xpos,(fail_split(1,ii)+fail_split(2,ii)),ratio_error(ii), 'k', 'LineStyle','none')
temp_text = {'n = '; [num2str(round(number_prch(ii)/1000)) 'k']};
text(xpos-width/2, (fail_split(1,ii)+fail_split(2,ii))+3, temp_text);
end
xlim([0.5 7*gap+3])
ylim([0 20])
ylabel('% of treatment failures ')
set(gca,'XTick',gap:gap:gap*8);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(params.new_order));
xtickangle(45);

end
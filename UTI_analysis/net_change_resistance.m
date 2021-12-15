function []= net_change_resistance(UTI_cases,params)
% function takes UTI_case stuc and optional params and makes bar chart of
% the % of early recurrences which changed susceptibiltiy to specific antiboitic, for all untreated cases and 
% all cases treated the the specific antibtioicu

%for all the antibiotics
for drug = 1:params.number_drugs 
   
   % initially sensitive or resistant 
   all_sensitive =  ismember(UTI_cases.SMP_Res(:,drug),params.sensitive_group) & UTI_cases.hasdiag ;
   all_resistant =  ismember(UTI_cases.SMP_Res(:,drug),params.resistant_group)  & UTI_cases.hasdiag ;
   
   % susceptibilty of recurrent infection
   all_sensitive_nextres =  UTI_cases.next_res(:,drug) == 1 | UTI_cases.SMP_Res(:,drug) == 2 ;
   all_resistant_nextres =  UTI_cases.next_res(:,drug) == 3 ;
   all_nextres = ismember(UTI_cases.next_res(:,drug),[1 2 3]);
   
   %num treated with the antibitoic
   all_sensitive_purchased2 =  nnz((all_sensitive | all_resistant) & UTI_cases.PCR_sameday(:,drug));
   sensitive_purchased_fails2 = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure);
   sensitive_purchased_nextresfails2 = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_nextres);
   
   sensitive_purchased_fails_SS = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_sensitive_nextres);
   sensitive_purchased_fails_SR = nnz(all_sensitive & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_resistant_nextres);
   
   %rate of failure for suscptibitiy-matched 
   ratio2(drug) = sensitive_purchased_fails2/all_sensitive_purchased2*100;
   number_prch(drug) = all_sensitive_purchased2;
   ratio_error2(drug) = sqrt((ratio2(drug)/100).*(1-(ratio2(drug)/100))./(all_sensitive_purchased2))*100;
   %split rate of failure by mode
   fail_split2(1:2,drug) = [sensitive_purchased_fails_SS/sensitive_purchased_nextresfails2 ; sensitive_purchased_fails_SR/sensitive_purchased_nextresfails2]*ratio2(drug);
   SR_error2(drug) = sqrt((fail_split2(2,drug)/100).*(1-(fail_split2(2,drug)/100))./(all_sensitive_purchased2))*100;
   
   %rate of failure for suscptibitiy-mismatched
   all_resistant_purchased2 = nnz((all_sensitive | all_resistant) & UTI_cases.PCR_sameday(:,drug));
   resistant_purchased_fails2 = nnz(all_resistant & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure);
   resistant_purchased_nextresfails2 = nnz(all_resistant & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_nextres);
   resistant_purchased_fails_RS = nnz(all_resistant & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_sensitive_nextres);
   resistant_purchased_fails_RR = nnz(all_resistant & UTI_cases.PCR_sameday(:,drug) & UTI_cases.treatfailure & all_resistant_nextres);
   
   ratio_res(drug) = resistant_purchased_fails2/all_resistant_purchased2*100;
   number_prch_res(drug) = all_resistant_purchased2;
   ratio_error_res(drug) = sqrt((ratio_res(drug)/100).*(1-(ratio_res(drug)/100))./(all_resistant_purchased2))*100;
   fail_split_res2(1:2,drug) = [resistant_purchased_fails_RR/resistant_purchased_nextresfails2 ; resistant_purchased_fails_RS/resistant_purchased_nextresfails2]*ratio_res(drug);
   SR_error_res2(drug) = sqrt((fail_split_res2(2,drug)/100).*(1-(fail_split_res2(2,drug)/100))./(all_resistant_purchased2))*100;
   
   %same for untreated
   all_sensitive_noprch2 =  nnz((all_sensitive | all_resistant) & UTI_cases.PCR_sameday(:,10));
   sensitive_noprch_fails2 = nnz(all_sensitive & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure );
   sensitive_noprch_nextresfails2 = nnz(all_sensitive & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure & all_nextres);
   
   sensitive_noprch_fails_SS = nnz(all_sensitive & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure & all_sensitive_nextres);
   sensitive_noprch_fails_SR = nnz(all_sensitive & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure & all_resistant_nextres);
   
   ratio_noprch2(drug) = sensitive_noprch_fails2/all_sensitive_noprch2*100;
   number_noprch(drug) = all_sensitive_noprch2;
   ratio_noprch_error(drug) = sqrt((ratio_noprch2(drug)/100).*(1-(ratio_noprch2(drug)/100))./(all_sensitive_noprch2))*100;
   fail_split_noprch2(1:2,drug) = [sensitive_noprch_fails_SS/sensitive_noprch_nextresfails2 ; sensitive_noprch_fails_SR/sensitive_noprch_nextresfails2]*ratio_noprch2(drug);
   SR_noprch_error2(drug) = sqrt((fail_split_noprch2(2,drug)/100).*(1-(fail_split_noprch2(2,drug)/100))./(all_sensitive_noprch2))*100;
   
   
   %all_resistant_noprch =  nnz(all_resistant & UTI_cases.PCR_sameday(:,10));
   all_resistant_noprch2 =  nnz((all_sensitive | all_resistant) & UTI_cases.PCR_sameday(:,10));
   resistant_noprch_fails2 = nnz(all_resistant & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure);
   resistant_noprch_nextresfails2 = nnz(all_resistant & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure & all_nextres);
   resistant_noprch_fails_RS = nnz(all_resistant & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure & all_sensitive_nextres);
   resistant_noprch_fails_RR = nnz(all_resistant & UTI_cases.PCR_sameday(:,10) & UTI_cases.treatfailure & all_resistant_nextres);
   
   ratio_noprch_res2(drug) = resistant_noprch_fails2/all_resistant_noprch2*100;
   number_noprch_res(drug) = all_resistant_noprch2;
   ratio_error_res(drug) = sqrt((ratio_noprch_res2(drug)/100).*(1-(ratio_noprch_res2(drug)/100))./(all_resistant_noprch2))*100;
   fail_split_noprch_res2(1:2,drug) = [resistant_noprch_fails_RR/resistant_noprch_nextresfails2 ; resistant_noprch_fails_RS/resistant_noprch_nextresfails2]*ratio_noprch_res2(drug);
   SR_error_noprch_res2(drug) = sqrt((fail_split_noprch_res2(2,drug)/100).*(1-(fail_split_noprch_res2(2,drug)/100))./(all_resistant_noprch2))*100;
     
   clear  all_sensitive all_sensitive_purchased sensitive_purchased_fails
end

%% plot the figure

width = 0.5;
gap = 1.8;
plot([0 7*gap+3], [0 0],'k')

for jj = 1:params.number_drugs
    ii = params.new_order(jj);
hold on
xpos = jj*gap-0.32;
hd3(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 fail_split_noprch2(2,ii) fail_split_noprch2(2,ii)],'w','EdgeColor','k', 'FaceColor', params.SR_color); % using MATLAB "patch" to establish the border
errorbar(xpos,fail_split_noprch2(2,ii)',SR_noprch_error2(ii), 'k', 'LineStyle','none')
%hatch(hda(ii),[45,4,1],'w')
hd4(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 -fail_split_noprch_res2(2,ii) -fail_split_noprch_res2(2,ii)],'w','EdgeColor','k', 'FaceColor', params.RS_color); % using MATLAB "patch" to establish the border
errorbar(xpos,-fail_split_noprch_res2(2,ii)',SR_error_noprch_res2(ii), 'k', 'LineStyle','none')
end

for jj = 1:params.number_drugs
    ii = params.new_order(jj);
hold on
xpos = jj*gap+0.32;
hd1(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 fail_split2(2,ii) fail_split2(2,ii)],'w','EdgeColor','k', 'FaceColor', params.SR_color); % using MATLAB "patch" to establish the border
errorbar(xpos,fail_split2(2,ii)',SR_error2(ii), 'k', 'LineStyle','none')
hd2(jj) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 -fail_split_res2(2,ii) -fail_split_res2(2,ii)],'w','EdgeColor','k', 'FaceColor', params.RS_color); % using MATLAB "patch" to establish the border
errorbar(xpos,-fail_split_res2(2,ii)',SR_error_res2(ii), 'k', 'LineStyle','none');
hatch(hd3(jj),[45,2.5,0.8],'k')
hatch(hd4(jj),[45,2.5,0.8],'k')
end
xlim([0.5 7*gap+3])
ylim([-1.5 7.5])
ylabel('% of patients who have an infection which changes resistance')

set(gca,'XTick',gap:gap:gap*8);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(params.new_order));
xtickangle(45);

end

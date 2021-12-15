function [] = initally_ecoli_gained_res_changed_bac(UTI_cases,params)
%all recurrent cases of initially E coli
% function takes UTI_case structure and optional params and calcuates the
% rate of change of species for initially E coli susceptibiltiy-matched treated cases which result in
% early recurrence, either S-S, or S-R.  

%%
% use periods for which relevant drug was routeenly measured
dates_to_use_start([1:4 6 8]) = min(UTI_cases.SamplingDate);
dates_to_use_start(5) = min(UTI_cases.SamplingDate)+7*321;
dates_to_use_start(7) = min(UTI_cases.SamplingDate)+7*293;
dates_to_use_end(1:7) = max(UTI_cases.SamplingDate);
dates_to_use_end(8) = min(UTI_cases.SamplingDate)+7*293;

% for gained resistance cases
%for all antibiotics
for drug = 1:params.number_drugs 
   dates_index =  find(UTI_cases.SamplingDate >= dates_to_use_start(drug) & UTI_cases.SamplingDate <= dates_to_use_end(drug));
   all_sensitive = (UTI_cases.SMP_Res(dates_index,drug) == 1 | UTI_cases.SMP_Res(dates_index,drug) == 2) & UTI_cases.hasdiag(dates_index) & UTI_cases.bug(dates_index) == 179 & UTI_cases.nBac(dates_index) == 1; % code 179 = E. coli
   all_resistant_nextres =  ismember(UTI_cases.next_res(dates_index,drug), params.resistant_group) ;
   sensitive_purchased_fails = nnz(all_sensitive & UTI_cases.PCR_sameday(dates_index,drug) & UTI_cases.treatfailure(dates_index) & all_resistant_nextres);
   
   % for the most common species
   sensitive_purchased_fails_stayedEcoli = nnz(all_sensitive & UTI_cases.PCR_sameday(dates_index,drug) & UTI_cases.treatfailure(dates_index) & all_resistant_nextres & UTI_cases.changed_bac(dates_index)==0);
   sensitive_purchased_fails_changed_kleb = nnz(all_sensitive & UTI_cases.PCR_sameday(dates_index,drug) & UTI_cases.treatfailure(dates_index) & all_resistant_nextres & UTI_cases.changed_bac(dates_index)==1 & UTI_cases.new_bac(dates_index)==225);
   sensitive_purchased_fails_changed_prot = nnz(all_sensitive & UTI_cases.PCR_sameday(dates_index,drug) & UTI_cases.treatfailure(dates_index) & all_resistant_nextres & UTI_cases.changed_bac(dates_index)==1 & UTI_cases.new_bac(dates_index)==280);
   sensitive_purchased_fails_changed_other = nnz(all_sensitive & UTI_cases.PCR_sameday(dates_index,drug) & UTI_cases.treatfailure(dates_index) & all_resistant_nextres & UTI_cases.changed_bac(dates_index)==1 & UTI_cases.new_bac(dates_index)~=280 & UTI_cases.new_bac(dates_index)~=225 & floor(UTI_cases.new_bac(dates_index))~=179 & UTI_cases.new_bac(dates_index)>0);
   
   ratio_stayed(drug) = sensitive_purchased_fails_stayedEcoli/sensitive_purchased_fails*100;
   ratio_kleb(drug) = sensitive_purchased_fails_changed_kleb/sensitive_purchased_fails*100;
   ratio_prot(drug) = sensitive_purchased_fails_changed_prot/sensitive_purchased_fails*100;
   ratio_other(drug) = sensitive_purchased_fails_changed_other/sensitive_purchased_fails*100;
   ratio_changed(drug) =  ratio_kleb(drug) + ratio_prot(drug) + ratio_other(drug);
   ratio_error(drug) = sqrt((ratio_changed(drug)/100).*(1-(ratio_changed(drug)/100))./(sensitive_purchased_fails))*100;
end

%%
total = sum([ratio_stayed; ratio_kleb; ratio_prot; ratio_other]);
ratio_stayed = ratio_stayed./total*100;
ratio_kleb = ratio_kleb./total*100;
ratio_prot = ratio_prot./total*100;
ratio_other = ratio_other./total*100;
gainedR = [ratio_other; ratio_prot; ratio_kleb];

%% remained sensitve 
for drug = 1:params.number_drugs 
   dates_index =  find(UTI_cases.SamplingDate >= dates_to_use_start(drug) & UTI_cases.SamplingDate <= dates_to_use_end(drug));
   all_sensitive = (UTI_cases.SMP_Res(dates_index,drug) == 1 | UTI_cases.SMP_Res(dates_index,drug) == 2) & UTI_cases.hasdiag(dates_index) & UTI_cases.bug(dates_index) == 179 & UTI_cases.nBac(dates_index) == 1;
   all_sensitive_nextres =  ismember(UTI_cases.next_res(dates_index,drug), params.sensitive_group) ;
   sensitive_purchased_fails = nnz(all_sensitive & UTI_cases.PCR_sameday(dates_index,drug) & UTI_cases.treatfailure(dates_index) & all_sensitive_nextres);
   sensitive_purchased_fails_stayedEcoli = nnz(all_sensitive & UTI_cases.PCR_sameday(dates_index,drug) & UTI_cases.treatfailure(dates_index) & all_sensitive_nextres & UTI_cases.changed_bac(dates_index)==0);
   sensitive_purchased_fails_changed_kleb = nnz(all_sensitive & UTI_cases.PCR_sameday(dates_index,drug) & UTI_cases.treatfailure(dates_index) & all_sensitive_nextres & UTI_cases.changed_bac(dates_index)==1 & UTI_cases.new_bac(dates_index)==225);
   sensitive_purchased_fails_changed_prot = nnz(all_sensitive & UTI_cases.PCR_sameday(dates_index,drug) & UTI_cases.treatfailure(dates_index) & all_sensitive_nextres & UTI_cases.changed_bac(dates_index)==1 & UTI_cases.new_bac(dates_index)==280);
   sensitive_purchased_fails_changed_other = nnz(all_sensitive & UTI_cases.PCR_sameday(dates_index,drug) & UTI_cases.treatfailure(dates_index) & all_sensitive_nextres & UTI_cases.changed_bac(dates_index)==1 & UTI_cases.new_bac(dates_index)~=280 & UTI_cases.new_bac(dates_index)~=225 & floor(UTI_cases.new_bac(dates_index))~=179 & UTI_cases.new_bac(dates_index)>0);
   
   ratio_stayed(drug) = sensitive_purchased_fails_stayedEcoli/sensitive_purchased_fails*100;
   ratio_kleb(drug) = sensitive_purchased_fails_changed_kleb/sensitive_purchased_fails*100;
   ratio_prot(drug) = sensitive_purchased_fails_changed_prot/sensitive_purchased_fails*100;
   ratio_other(drug) = sensitive_purchased_fails_changed_other/sensitive_purchased_fails*100;
   ratio_changed(drug) =  ratio_kleb(drug) + ratio_prot(drug) + ratio_other(drug);
   ratio_error2(drug) = sqrt((ratio_changed(drug)/100).*(1-(ratio_changed(drug)/100))./(sensitive_purchased_fails))*100;
   
end
%%
total = sum([ratio_stayed; ratio_kleb; ratio_prot; ratio_other]);
ratio_stayed = ratio_stayed./total*100;
ratio_kleb = ratio_kleb./total*100;
ratio_prot = ratio_prot./total*100;
ratio_other = ratio_other./total*100;
remainedS = [ratio_other; ratio_prot; ratio_kleb];

%% plot the fig
%reorder for fig
gainedR = gainedR(:,params.new_order);
remainedS =remainedS(:,params.new_order);

width = 0.5;
gap = 1.8;
plot([0 7*gap+3], [0 0],'k')
new_kleb_col = [0.5, 0.5, 0.5];
new_prot_col = [0.8, 0.8, 0.8];
new_other_col = [0.2, 0.2, 0.2];
for ii = 1:8
hold on
xpos = ii*gap-0.32;
hd1(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 remainedS(1,ii) remainedS(1,ii)],'w','EdgeColor','k', 'FaceColor', new_other_col); % using MATLAB "patch" to establish the border
hd2(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[remainedS(1,ii)  remainedS(1,ii) (remainedS(1,ii)+remainedS(2,ii)) (remainedS(1,ii)+remainedS(2,ii))],'w','EdgeColor','k', 'FaceColor', new_prot_col); % using MATLAB "patch" to establish the border
hd3(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[sum(remainedS(1:2,ii)) sum(remainedS(1:2,ii)) sum(remainedS(1:3,ii)) sum(remainedS(1:3,ii))],'w','EdgeColor','k', 'FaceColor', new_kleb_col ); % using MATLAB "patch" to establish the border
errorbar(xpos,sum(remainedS(1:3,ii))',ratio_error2(ii), 'k', 'LineStyle','none')

end
for ii = 1:8
hold on
xpos = ii*gap+0.32;
hd1(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 gainedR(1,ii) gainedR(1,ii)],'w','EdgeColor','k', 'FaceColor', new_other_col); % using MATLAB "patch" to establish the border
hd2(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[gainedR(1,ii)  gainedR(1,ii) (gainedR(1,ii)+gainedR(2,ii)) (gainedR(1,ii)+gainedR(2,ii))],'w','EdgeColor','k', 'FaceColor', new_prot_col); % using MATLAB "patch" to establish the border
hd3(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[sum(gainedR(1:2,ii)) sum(gainedR(1:2,ii)) sum(gainedR(1:3,ii)) sum(gainedR(1:3,ii))],'w','EdgeColor','k', 'FaceColor', new_kleb_col ); % using MATLAB "patch" to establish the border
errorbar(xpos,sum(gainedR(1:3,ii))',ratio_error(ii), 'k', 'LineStyle','none')

end

xlim([0.5 7*gap+3])
ylim([0 80])
% xlim([0 8.5])
% ylim([0 80])
ylabel({'% of initially {\it E. coli} infections which'; ...
    'gained resistance that were caused by';...
    'reinfection different species'})
xlabel('treatment antibiotic')
set(gca,'XTick',width:gap:gap*8);
set(gca,'xticklabel', UTI_cases.SMP_Res_drug_names([1 2 8 3:7]));
xtickangle(45);
legend('Other species','{\it P. mirabilis}','{\it K. pneumoniae}')

function [] = mode_of_reccurence_pie_charts_wounds(wound_cases, params)
% function takes wound_cases strcut and optional params and makes pie charts of 
% the % of early recurrences and their mode of recurrence for all antibiotic 
% treated cases.


%define fig properties
figure
Names = wound_cases.SMP_Res_drug_names;
set(gcf,'color','w', 'name','Fig. 1H');
left_val = 5;
bottom_val = 5;
width_val = 17;
height_val = 10;
set(gcf,'units','centimeters','Position',[left_val bottom_val width_val height_val]);


for drug = 1:params.number_drugs 
   total_prch1(drug) = nnz(wound_cases.PCR_sameday(:,drug));
   num_fails(drug) = nnz(wound_cases.PCR_sameday(:,drug)  & wound_cases.treatfailure); 

   total_prch_matched(drug) = nnz(wound_cases.PCR_sameday(:,drug) & ismember(wound_cases.SMP_Res(:,drug),[1 2]));
   total_prch_mismatched(drug) = nnz(wound_cases.PCR_sameday(:,drug)  & ismember(wound_cases.SMP_Res(:,drug),3));
 
   total_prch_matched_fails(drug) = nnz(wound_cases.PCR_sameday(:,drug) & ismember(wound_cases.SMP_Res(:,drug),[1 2]) & wound_cases.treatfailure);
   total_prch_mismatched_fails(drug) = nnz(wound_cases.PCR_sameday(:,drug)  & ismember(wound_cases.SMP_Res(:,drug),3) & wound_cases.treatfailure);
    
   total_prch_matched_clear(drug) = nnz(wound_cases.PCR_sameday(:,drug)  & ismember(wound_cases.SMP_Res(:,drug),[1 2]) & ~wound_cases.treatfailure);
   total_prch_mismatched_clear(drug) = nnz(wound_cases.PCR_sameday(:,drug) & ismember(wound_cases.SMP_Res(:,drug),3) & ~wound_cases.treatfailure);
end

num_total_prch = sum(total_prch_matched) + sum(total_prch_mismatched);
num_early_rec = sum(num_fails)/sum(total_prch1);
fprintf('The total number of treated infections is %i\n',num_total_prch);
fprintf('The percentage number of early recurrences is %2.1f\n',num_early_rec*100);

num_early_rec = (sum(total_prch_mismatched_fails)+sum(total_prch_matched_fails));
num_prch_matched_clear = sum(total_prch_matched_clear);
num_prch_mismatched_clear = sum(total_prch_mismatched_clear);

per_early_rec = num_early_rec/(sum(total_prch_mismatched)+sum(total_prch_matched));
per_prch_matched_clear = num_prch_matched_clear/(sum(total_prch_mismatched)+sum(total_prch_matched));
per_prch_mismatched_clear = num_prch_mismatched_clear/(sum(total_prch_mismatched)+sum(total_prch_matched));

t = tiledlayout(3,5, 'Padding','tight');
nexttile(1,[2 2]);
pie([per_early_rec per_prch_matched_clear per_prch_mismatched_clear ], ...
    {sprintf('early reccurrence \n %d',num_early_rec),...
    sprintf('S->0 \n %d',num_prch_matched_clear),...
    sprintf('R->0 \n %d',num_prch_mismatched_clear)});
ax = gca();
ax.Colormap = [[0.32 0.32 0.32];  [0.74, 0.74, 0.98]; [0.98, 0.79, 0.79]];
title('all treated UTIs')


%
no_measurement = wound_cases.SMP_Res == 0;
no_next_measurement = wound_cases.next_res == 0 | isnan(wound_cases.next_res);
no_measurement_either_all = no_measurement | no_next_measurement;
no_measurement_either = no_measurement_either_all(wound_cases.treatfailure == 1,:);
CurrentRes = wound_cases.SMP_Res(wound_cases.treatfailure == 1,:) ;
CurrentRes = ismember(CurrentRes,params.resistant_group);
NextRes = wound_cases.next_res(wound_cases.treatfailure == 1,:);
NextRes = ismember(NextRes,params.resistant_group);
SS = ~CurrentRes & ~NextRes & ~no_measurement_either;
SR = ~CurrentRes & NextRes & ~no_measurement_either;
RS = CurrentRes & ~NextRes & ~no_measurement_either;
RR = CurrentRes & NextRes & ~no_measurement_either;
drug_prch_fails = wound_cases.PCR_sameday(wound_cases.treatfailure,:);

for ii = 1:params.number_drugs
SS_drugprch = SS(drug_prch_fails(:,ii),ii);
SR_drugprch = SR(drug_prch_fails(:,ii),ii);
RS_drugprch = RS(drug_prch_fails(:,ii),ii);
RR_drugprch = RR(drug_prch_fails(:,ii),ii);

SSnum(ii) = length(SS_drugprch(SS_drugprch));
SRnum(ii) = length(SR_drugprch(SR_drugprch));
RSnum(ii) = length(RS_drugprch(RS_drugprch));
RRnum(ii) = length(RR_drugprch(RR_drugprch));
total(ii) = SSnum(ii) + SRnum(ii) + RSnum(ii) + RRnum(ii) ;
fraction = [ RRnum(ii)  RSnum(ii) SRnum(ii) SSnum(ii)]./total(ii);

ax = nexttile(10+ii);
pie(fraction);
ax.Colormap = [params.RR_color;  params.RS_color; params.SR_color;  params.SS_color];
TextChildren =  findobj(ax,'Type','text');
set(TextChildren,'visible','off')

title(Names{ii})
savefraction(:,ii) = fraction;

clear fraction SS_drugprch SR_drugprch RS_drugprch RR_drugprch
end

fraction = [sum(RRnum) sum(RSnum) sum(SRnum) sum(SSnum) ]./(sum(SSnum) + sum(SRnum) + sum(RSnum) + sum(RRnum));
nexttile(3,[2 2])
pie(fraction);
ax = gca();
ax.Colormap = [params.RR_color;  params.RS_color; params.SR_color;  params.SS_color];

end

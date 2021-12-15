%% perform UTI analysis and makes figures for 
% Minimizing antibiotic-induced emergence of antibiotic resistance in
% bacterial infections
% Mathew Stracy 2021
%% loads the UTI_cases structure. Each line in the struct represents a UTI case
clearvars -except UTI_cases
if ~exist('UTI_cases') % load the data
    load('UTI_cases.mat');
end
if ~isfolder('Tables')
      mkdir('Tables')
end

%% parameters
params.number_drugs = 8; 
% sensitivity grouping Sensitive = 1; Intermediate = 2; Resistant = 3;
params.sensitive_group = [1 2];
SIR= 1:3; params.resistant_group = SIR(~ismember(SIR, params.sensitive_group));
clear SIR
% figure colors
params.SS_color = [0.1 0.13 0.5]; params.SR_color = [0 0.95 0.95];
params.RR_color = [0.5 0 0.5]; params.RS_color = [0.95 0 0.95];
params.new_order = [1 2 8 3:7]; %drug order for figures

%% generate Table S1 patient demographics 
Patient_demographics_table_agecat(UTI_cases);

%% Fig 1B
mismatched_vs_matched_recurrence_rate(UTI_cases,params)

%% Fig 1C
mode_of_reccurence_pie_charts(UTI_cases,params)

%% Fig 1D 
figure
set(gcf,'color','w','name','Fig. 1D','units','centimeters','Position',[1 1 20 7])
risk_reccurence_correctly_prescribed_treat_v_untreat(UTI_cases,params)

%% Fig S3
figure
set(gcf,'color','w', 'name','Fig. S3','units','centimeters','Position',[1 1 17 7.5]);
net_change_resistance(UTI_cases,params)

%% Fig 1E rate of recuurence by day 
figure
set(gcf,'color','w','name','Fig. 1E','units','centimeters','Position',[1 1 11 6]);
risk_recurrence_by_day(UTI_cases,params)
drawnow

%% Fig S6 Mode of recurrence by day for treated and untreated UTIs. Mode of recurrence by day for treated and untreated UTIs. 
figure
set(gcf,'color','w', 'name','Fig. S6 Mode of recurrence','units','centimeters','Position',[1 1 16 26]);
mode_of_recurrence_treated_untreated(UTI_cases,params)
drawnow

%% Fig 1F matrix UTIs
figure
set(gcf,'color','w', 'name','Fig. 1F Susceptibiltiy change matrix', 'units','centimeters','Position',[1 1 12 10]);
susceptibiltiy_change_matrix2(UTI_cases,params)
drawnow

%% Fig 2E
figure;
set(gcf,'color','w', 'name','Fig. 2F', 'units','centimeters','Position',[1 1 18 10]);
initally_ecoli_gained_res_changed_bac(UTI_cases,params)

%% Fig 2G
figure
set(gcf,'color','w', 'name','Fig. 2G UTIs', 'units','centimeters','Position',[1 1 5 9]);
changed_res_bac_alldrugs(UTI_cases,params)

%% Fig 3B
figure;
set(gcf,'color','w', 'name','Fig. 3B Adjusted OR');
correct_treatment_adjusted_OR(UTI_cases,params)
drawnow

%% Fig 3DEF
correct_treatment_fails_reccomend_drug(UTI_cases,params)
drawnow
%% Supp Fig S5 
figure;
set(gcf,'color','w')
clf; set(gcf,'name','Fig. S5 Adjusted risk of recurrence and emergence of resistance with and without antibiotic treatment');
set(gcf,'units','centimeters','Position',[1 1 10 30]);
Adjusted_risk_of_recurrence_treat_notreat(UTI_cases,params)

%% Supp Fig S8
figure;
set(gcf,'color','w', 'name','Fig. S8 Rates of resistance by species UTIs',...
'units','centimeters','Position',[1 1 25 18]);
rate_res_by_species(UTI_cases,params);

%% Supp Fig S10
Prev_res_OR_byYear2(UTI_cases,params)

%% Supp Fig S11
figure;
set(gcf,'color','w', 'name','Fig. S11', 'units','centimeters','Position',[1 1 25 18]);
Adjucted_OR_PURCHandRES(UTI_cases,params)
Prev_res_OR_exactly1(UTI_cases,params)

%% Supp Fig S14BC
correct_treatment_fails_reccomend_drug_v9_2testperiods(UTI_cases,params)
drawnow
%% Supp Fig 9
species_change_matrix_gained_res(UTI_cases,params)

%% Supp Fig S16 intermediate grouped with resistant
figure;
set(gcf,'color','w', 'name','Fig. S16', 'units','centimeters','Position',[1 1 25 10]);
subplot(1,2,1)
correct_treatment_adjusted_OR_S_IR(UTI_cases,params);
correct_treatment_fails_reccomend_drug_S_IR(UTI_cases,params)
drawnow

%% Supp Table S4
species_prevalence_table(UTI_cases)

%% save the figures
if ~isfolder('Figures')
      mkdir('Figures')
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'name');
  saveas(FigHandle,['Figures/' FigName '.fig'])
end
%% perform wound infection analysis and makes figures for 
% Minimizing antibiotic-induced emergence of antibiotic resistance in
% bacterial infections
% Mathew Stracy 2021

%% loads the wound_cases structure. Each line in the struct represents a wound infection case
clearvars -except wound_cases
if ~exist('wound_cases') % load the data
    load('wound_cases.mat');
end
if ~isfolder('Tables')
      mkdir('Tables')
end

%% parameters
params.number_drugs = 5; 
% sensitivity grouping Sensitive = 1; Intermediate = 2; Resistant = 3;
params.sensitive_group = [1 2];
SIR= 1:3; params.resistant_group = SIR(~ismember(SIR, params.sensitive_group));
clear SIR
% figure colors
params.SS_color = [ 25, 32, 128]/255; params.SR_color = [ 0 0.95 0.95];
params.RR_color = [ 128 0 128]/255; params.RS_color = [ 0.95 0 0.95];

%% generate Table S1 patient demographics 
Patient_demographics_table_wounds2(wound_cases,params);

%% Fig 1G
figure
set(gcf,'color','w','name','Fig. 1G', 'units','centimeters','Position',[2 2 7 10]);
mismatched_vs_mathced_recurrence_rate_wounds(wound_cases,params)

%% Fig 1H
mode_of_reccurence_pie_charts_wounds(wound_cases,params)

%% Table S4
species_prevalence_table(wound_cases)

%% Fig 1H
figure
set(gcf,'color','w', 'name','Fig. 1H', 'units','centimeters','Position',[1 1 10 10]);
susceptibiltiy_change_matrix_wounds(wound_cases,params)

%% Fig S4
figure
set(gcf,'color','w', 'name','Fig. S4', 'units','centimeters','Position',[1 1 12 9]);
risk_recurrence_by_day_wounds(wound_cases,params)

%% Fig 2G
figure
set(gcf,'color','w', 'name','Fig. 2G wounds', 'units','centimeters','Position',[1 1 5 9]);
changed_res_bac_alldrugs_wounds(wound_cases,params)

%% Fig S8 wounds
figure;
set(gcf,'color','w', 'name','Fig. S8B Rates of resistance by species wounds',...
'units','centimeters','Position',[1 1 25 18]);
rate_res_by_species_wounds(wound_cases, params);

%% Fig 3C
figure;
set(gcf,'color','w', 'name','Fig. 4C',...
'units','centimeters','Position',[1 1 25 18]);
OR_wounds(wound_cases, params);

%% Fig 3G
set(gcf,'color','w', 'name','Fig. S12 wounds',...
'units','centimeters','Position',[1 1 25 18]);
reccomend_drugs(wound_cases, params)

%%
species_change_matrix_gained_res_wound(wound_cases,params)
%%
if ~isfolder('Figures')
      mkdir('Figures')
end
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'name');
  saveas(FigHandle,['Figures/' FigName '.fig'])
end
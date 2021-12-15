function [] = species_prevalence_table(UTI_cases)
%function makes table of the species prevelance for all UTI cases
A = cell2mat(UTI_cases.bug_all);
[counts,bugs] = groupcounts(A);
species_prevalence_table = table(bugs,counts);
total = sum(counts);
for ii = 1:height(species_prevalence_table)
   names(ii,1) = UTI_cases.Bugs.Name( UTI_cases.Bugs.Code == species_prevalence_table.bugs(ii));
   percent(ii,1) = counts(ii)/total*100;
end
species_prevalence_table.Names = names;
species_prevalence_table.percent = round(percent);
species_prevalence_table_wounds = sortrows(species_prevalence_table,'counts','descend');
species_prevalence_table_wounds = species_prevalence_table_wounds(1:10,:);
filename = 'Tables/TableS4_species_prevalence_table_UTIs.xlsx';
writetable(species_prevalence_table_wounds,filename);
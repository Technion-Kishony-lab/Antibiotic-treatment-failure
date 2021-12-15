function [] = changed_res_bac_alldrugs_wounds(wound_cases,params)
%function takes the wound_cases structure and optional params and calulates
% the rate of change of species in gained-resistance early recurrences 

treatfailure = wound_cases.treatfailure;
clear ratio_error
for drug = 1:params.number_drugs 
    
%next res is resistant
index = cellfun(@(x) x==3, wound_cases.new_RES(:,drug),'UniformOutput',false);
%only resistant next wound_cases.Bugs
resistant_nextres_wound_cases.Bugs = cellfun(@(v,x)v(x),wound_cases.new_bug,index,'UniformOutput',false);
%are any of the next resistant wound_cases.Bugs not those in original infection
[M2,X2] = cellfun(@ismember,resistant_nextres_wound_cases.Bugs,wound_cases.bug_all,'UniformOutput',false);

resistant_nextres_wound_cases.Bugs_members = cellfun(@(v,x)v(x),wound_cases.new_bug,M2,'UniformOutput',false);
M3 = cellfun(@not,M2,'UniformOutput',false);
resistant_nextres_wound_cases.Bugs_notmembers = cellfun(@(v,x)v(x),wound_cases.new_bug,M3,'UniformOutput',false);

initially_sensitive = (wound_cases.SMP_Res(:,drug) == 1 | wound_cases.SMP_Res(:,drug) == 2);

to_use = find(initially_sensitive & treatfailure & wound_cases.PCR_sameday(:,drug));

changed_species = cell2mat( resistant_nextres_wound_cases.Bugs_notmembers(to_use));
[GC,GR] = groupcounts(changed_species) ;
table(wound_cases.Bugs.Name(ismember(wound_cases.Bugs.Code, GR)),GC);

same_species = cell2mat( resistant_nextres_wound_cases.Bugs_members(to_use));
[GC,GR] = groupcounts(same_species) ;
table(wound_cases.Bugs.Name(ismember(wound_cases.Bugs.Code, GR)),GC);

changed(drug) = length(changed_species);
stayed(drug) = length(same_species);

ratio_changed(drug) = changed(drug)./(changed(drug) + stayed(drug));
ratio_error(drug) = sqrt((ratio_changed(drug)).*(1-(ratio_changed(drug)))./(changed(drug)+stayed(drug)));
 nums(drug) = (changed(drug) + stayed(drug))  ;

end

ratio_changed2 = sum(changed)./(sum(changed) + sum(stayed));
ratio_error2 = sqrt((ratio_changed2).*(1-(ratio_changed2))./(sum(changed) + sum(stayed)));
% ratio_changed2(2) = 0.7466;
% ratio_error2(2) = 0.0253;
bar(ratio_changed2, 'FaceColor',params.SR_color)
xticklabels({'Wounds'})
xtickangle(45)
hold on
errorbar(1:length(ratio_changed2),ratio_changed2,ratio_error2, 'k', 'LineStyle','none')
ylabel('% of SR caused by different species')
ylim([0 1])


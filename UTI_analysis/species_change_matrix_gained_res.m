function [] = species_change_matrix_gained_res(UTI_cases,params)
%function takes the UTI_case structure and optional params and calulates
% matrix of the rate of change of species in gained-resistance early recurrences 

%% matrix UTI gained res

bugs_order = [1 381 292 110 280 225 179]; % most common UTI species, 1 = other
clear names
for ii = 1:length(bugs_order)
   names(ii) = UTI_cases.Bugs.Name(find(UTI_cases.Bugs.Code == bugs_order(ii))); 
end
names(1) = {'other'};
not_members = 1:1000;
not_members(bugs_order) = [];
clear mat_changed
treatfailure = UTI_cases.treatfailure;% 
mat_changed  = zeros(length(bugs_order),length(bugs_order),params.number_drugs );

for drug = 1:params.number_drugs 
    
%next res is resistant
index = cellfun(@(x) ismember(x,params.resistant_group), UTI_cases.new_RES(:,drug),'UniformOutput',false);
%only resistant next bugs
resistant_nextres_bugs = cellfun(@(v,x)v(x),UTI_cases.new_bug,index,'UniformOutput',false);

initially_sensitive = (UTI_cases.SMP_Res(:,drug) == 1 | UTI_cases.SMP_Res(:,drug) == 2);
to_use = find(initially_sensitive & treatfailure & UTI_cases.PCR_sameday(:,drug));

for inital_bug = 1:length(bugs_order)
    if inital_bug == 1 %all other species
    initial_true =  cellfun(@(x) ismember(x,not_members),UTI_cases.bug_all,'UniformOutput',false);
    initial_true = cellfun(@sum,initial_true)>0;    
    else
    initial_true = (cellfun(@(x) x==bugs_order(inital_bug),UTI_cases.bug_all,'UniformOutput',false));
    initial_true = cellfun(@sum,initial_true)>0;
    end
    for new_bug = 1:length(bugs_order)
         if new_bug == 1
    recurrent_true =  cellfun(@(x) ismember(x,not_members),resistant_nextres_bugs,'UniformOutput',false);
    recurrent_true = cellfun(@sum,recurrent_true)>0;    
    else
    recurrent_true = (cellfun(@(x) x==bugs_order(new_bug),resistant_nextres_bugs,'UniformOutput',false));
    recurrent_true = cellfun(@sum,recurrent_true)>0; 
         end
    mat_changed(inital_bug,new_bug,drug) = nnz(initial_true(to_use) & recurrent_true(to_use));  
    end
end

end

%%
figure;
set(gcf,'color','w', 'name','Fig. S9 UTIs drugs', 'units','centimeters','Position',[1 1 18 10]);
tiledlayout(2,4, 'TileSpacing','Compact', 'Padding', 'none');

for drug = 1:params.number_drugs
nexttile
temp = mat_changed(:,:,drug);
for ii = 1:length(bugs_order)
temp(ii,:) = temp(ii,:)./sum(temp(ii,:));
end
image(fliplr(temp*256))
title( UTI_cases.SMP_Res_drug_names(drug));
yticklabels([])
xticklabels([])
axis image

end

figure;
set(gcf,'color','w', 'name','Fig. S9 UTIs', 'units','centimeters','Position',[1 1 25 10]);
tiledlayout(1,2, 'TileSpacing','Compact', 'Padding', 'none');
nexttile
mat_changed_all = sum(mat_changed, 3);
temp = mat_changed_all;
for ii = 1:length(bugs_order)
temp(ii,:) = temp(ii,:)./sum(temp(ii,:));
end
image(fliplr(temp*256))
title('all antibiotics')
xticks(1:length(bugs_order))
yticklabels(names)
xtickangle(45)
yticks(1:length(bugs_order))
xticklabels(fliplr(names))
axis image
ylabel('Initial species')
xlabel('Gained resistance species')
%colorbar

%% matrix UTI remained sensitive
mat_changedSS  = zeros(length(bugs_order),length(bugs_order),params.number_drugs );

for drug = 1:params.number_drugs 
    
index = cellfun(@(x) ismember(x,params.sensitive_group), UTI_cases.new_RES(:,drug),'UniformOutput',false);
%only sensitive next bugs
resistant_nextres_bugs = cellfun(@(v,x)v(x),UTI_cases.new_bug,index,'UniformOutput',false);

initially_sensitive = (UTI_cases.SMP_Res(:,drug) == 1 | UTI_cases.SMP_Res(:,drug) == 2);
to_use = find(initially_sensitive & treatfailure & UTI_cases.PCR_sameday(:,drug));

for inital_bug = 1:length(bugs_order)
    if inital_bug == 1
    initial_true =  cellfun(@(x) ismember(x,not_members),UTI_cases.bug_all,'UniformOutput',false);
    initial_true = cellfun(@sum,initial_true)>0;    
    else
    initial_true = (cellfun(@(x) x==bugs_order(inital_bug),UTI_cases.bug_all,'UniformOutput',false));
    initial_true = cellfun(@sum,initial_true)>0;
    end
    for new_bug = 1:length(bugs_order)
         if new_bug == 1
    recurrent_true =  cellfun(@(x) ismember(x,not_members),resistant_nextres_bugs,'UniformOutput',false);
    recurrent_true = cellfun(@sum,recurrent_true)>0;    
    else
    recurrent_true = (cellfun(@(x) x==bugs_order(new_bug),resistant_nextres_bugs,'UniformOutput',false));
    recurrent_true = cellfun(@sum,recurrent_true)>0; 
         end
    mat_changedSS(inital_bug,new_bug,drug) = nnz(initial_true(to_use) & recurrent_true(to_use));  
    end
end


end

nexttile
mat_changed_all = sum(mat_changedSS, 3);
temp = mat_changed_all;
for ii = 1:length(bugs_order)
temp(ii,:) = temp(ii,:)./sum(temp(ii,:));
end
image(fliplr(temp*256))
title('all antibiotics')
xticks(1:length(bugs_order))
yticklabels(names)
xtickangle(45)
yticks(1:length(bugs_order))
xticklabels(fliplr(names))
axis image
ylabel('Initial species')
xlabel('Remained sensitive species')
%colorbar
end
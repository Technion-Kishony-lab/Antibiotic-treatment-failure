function [] = species_change_matrix_gained_res_wound(wound_cases,params)
%% matrix wound gained res

bugs_order = [1 5 223 110 242 166 225 280 179 292 347]; %wounds
clear names
for ii = 1:length(bugs_order)
   names(ii) = wound_cases.Bugs.Name(find(wound_cases.Bugs.Code == bugs_order(ii))); 
end
names(1) = {'other'};
not_members = 1:1000;
not_members(bugs_order) = [];
clear mat_changed
treatfailure = wound_cases.treatfailure;%   wound_cases.next_SMP_days<29 & wound_cases.next_SMP_days>4;
mat_changed  = zeros(length(bugs_order),length(bugs_order),params.number_drugs );

for drug = 1:params.number_drugs 
    
%next res is resistant
index = cellfun(@(x) ismember(x,params.resistant_group), wound_cases.new_RES(:,drug),'UniformOutput',false);

%only resistant next bugs
resistant_nextres_bugs = cellfun(@(v,x)v(x),wound_cases.new_bug,index,'UniformOutput',false);

initially_sensitive = (wound_cases.SMP_Res(:,drug) == 1 | wound_cases.SMP_Res(:,drug) == 2);
to_use = find(initially_sensitive & treatfailure & wound_cases.PCR_sameday(:,drug));

for inital_bug = 1:length(bugs_order)
    if inital_bug == 1
    initial_true =  cellfun(@(x) ismember(x,not_members),wound_cases.bug_all,'UniformOutput',false);
    initial_true = cellfun(@sum,initial_true)>0;    
    else
    initial_true = (cellfun(@(x) x==bugs_order(inital_bug),wound_cases.bug_all,'UniformOutput',false));
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

figure;
set(gcf,'color','w', 'name','Fig. S9 wounds', 'units','centimeters','Position',[1 1 25 10]);
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

%%
mat_changedSS  = zeros(length(bugs_order),length(bugs_order),params.number_drugs );

for drug = 1:params.number_drugs 
    
index = cellfun(@(x) ismember(x,params.sensitive_group), wound_cases.new_RES(:,drug),'UniformOutput',false);
%only sensitive next bugs
resistant_nextres_bugs = cellfun(@(v,x)v(x),wound_cases.new_bug,index,'UniformOutput',false);

initially_sensitive = (wound_cases.SMP_Res(:,drug) == 1 | wound_cases.SMP_Res(:,drug) == 2);
to_use = find(initially_sensitive & treatfailure & wound_cases.PCR_sameday(:,drug));

for inital_bug = 1:length(bugs_order)
    if inital_bug == 1
    initial_true =  cellfun(@(x) ismember(x,not_members),wound_cases.bug_all,'UniformOutput',false);
    initial_true = cellfun(@sum,initial_true)>0;    
    else
    initial_true = (cellfun(@(x) x==bugs_order(inital_bug),wound_cases.bug_all,'UniformOutput',false));
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
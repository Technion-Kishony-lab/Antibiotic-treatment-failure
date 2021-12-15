function [] = susceptibiltiy_change_matrix_wounds(wound_cases, params)
% function takes UTI_case structure and optional params and makes matric plot of 
% rate of change resistance for each antiboitic suceptibiltiy and each
% antiboitic treatment
 
%% net change in resistance to test antibiotic for each treatement antibiotic
% use periods for which relevant drug was routeenly measured
dates_to_use_start([1:4 6 8]) = min(wound_cases.SamplingDate);
dates_to_use_start(5) = min(wound_cases.SamplingDate)+7*321;
dates_to_use_start(7) = min(wound_cases.SamplingDate)+7*293;
dates_to_use_end(1:7) = max(wound_cases.SamplingDate);
dates_to_use_end(8) = min(wound_cases.SamplingDate)+7*293;


num_gained_resistance_to_test = zeros(params.number_drugs);
total_num_treated = zeros(params.number_drugs);
num_lost_resistance_to_test = zeros(params.number_drugs);


for drug = 1:params.number_drugs 
    total_prch(drug) = nnz(wound_cases.PCR_sameday(:,drug));
    dates_index_drug =  find(wound_cases.SamplingDate >= dates_to_use_start(drug) & wound_cases.SamplingDate <= dates_to_use_end(drug));       
    for drug_to_test = 1:params.number_drugs 
    dates_index_test =  find(wound_cases.SamplingDate >= dates_to_use_start(drug_to_test) & wound_cases.SamplingDate <= dates_to_use_end(drug_to_test));
    dates_index =  intersect(dates_index_test , dates_index_drug); 
    all_sensitive_test = (wound_cases.SMP_Res(dates_index,drug_to_test) == 1 | wound_cases.SMP_Res(dates_index,drug_to_test) == 2)  ;
    all_resistant_test = wound_cases.SMP_Res(dates_index,drug_to_test) == 3   ;
    all_currentres = ismember(wound_cases.SMP_Res(dates_index,drug),[1 2]) & ismember(wound_cases.SMP_Res(dates_index,drug_to_test),[1 2 3]);  
    all_sensitive_nexttestres =  wound_cases.next_res(dates_index,drug_to_test) == 1 | wound_cases.SMP_Res(dates_index,drug_to_test) == 2 ;
    all_resistant_nexttestres =  wound_cases.next_res(dates_index,drug_to_test) == 3 ;
    all_nextres = ismember(wound_cases.next_res(dates_index,drug),[1 2 3]) & ismember(wound_cases.next_res(dates_index,drug_to_test),[1 2 3]);
    gained_resistance_to_test = all_sensitive_test & all_resistant_nexttestres & wound_cases.PCR_sameday(dates_index,drug);
    lost_resistance_to_test = all_resistant_test & all_sensitive_nexttestres & wound_cases.PCR_sameday(dates_index,drug);   
    total_num_treated(drug,drug_to_test) = nnz(wound_cases.PCR_sameday(dates_index,drug) & all_currentres);
    num_gained_resistance_to_test(drug,drug_to_test) = nnz(gained_resistance_to_test & wound_cases.treatfailure(dates_index) & all_nextres & all_currentres);
    num_lost_resistance_to_test(drug,drug_to_test) = nnz(lost_resistance_to_test & wound_cases.treatfailure(dates_index) & all_nextres & all_currentres);
    end
end 



% net change in resistance
rate_gained_resistance_to_test = num_gained_resistance_to_test./total_num_treated*100;
rate_lost_resistance_to_test = num_lost_resistance_to_test./total_num_treated*100;
rate_gained_resistance_to_test(rate_gained_resistance_to_test == Inf) = 0;
rate_lost_resistance_to_test(rate_lost_resistance_to_test == Inf) = 0;
net_change_res = rate_gained_resistance_to_test-rate_lost_resistance_to_test;
net_change_res(isnan(rate_lost_resistance_to_test)) = 0;

%%
% for colormap
n = 20;
lost_res_num = 15;
gained_res_num = 55;
total_num = gained_res_num + lost_res_num;
cyan = [ linspace(0, 1, gained_res_num)' ones(gained_res_num,1) ones(gained_res_num,1)];
magenta = [  ones(gained_res_num,1) linspace(0, 1, gained_res_num)' ones(gained_res_num,1)];
magenta = flipud(magenta);
cyan_white_magenta = [cyan; magenta(1:lost_res_num,:)];
colormap(flipud(cyan_white_magenta));

image([rot90(net_change_res)]*n+lost_res_num )
xlabel('treatment drug')
ylabel('focal drug')
colormap(flipud(cyan_white_magenta))
xticks([1:1:params.number_drugs+1])
set(gca,'xticklabel',[wound_cases.SMP_Res_drug_names; 'Untreated'])
xtickangle(45);
yticks([1:1:params.number_drugs])
set(gca,'yticklabel',flipud(wound_cases.SMP_Res_drug_names))
cbh = colorbar ; %Create Colorbar
cbh.Ticks = linspace(0, length(colormap), 8) ; %Create 8 ticks from zero to 1
ticknums = (linspace(-lost_res_num, gained_res_num, 8))/n;
cbh.TickLabels = num2cell(ticknums) ; 
set(get(cbh,'label'),'string',{'net % of treated cases which changed';...
    'resistance to focal drug'});
axis image

end
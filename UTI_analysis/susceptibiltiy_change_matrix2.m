function [] = susceptibiltiy_change_matrix(UTI_cases,params)
% function takes UTI_case structure and optional params and makes matric plot of 
% rate of change resistance for each antiboitic suceptibiltiy and each
% antiboitic treatment
 
%% sensitive chance of failure
% use periods for which relevant drug was routeenly measured
dates_to_use_start([1:4 6 8]) = min(UTI_cases.SamplingDate);
dates_to_use_start(5) = min(UTI_cases.SamplingDate)+7*321; 
dates_to_use_start(7) = min(UTI_cases.SamplingDate)+7*293;
dates_to_use_end(1:7) = max(UTI_cases.SamplingDate);
dates_to_use_end(8) = min(UTI_cases.SamplingDate)+7*293; 

%preallocate
num_gained_resistance_to_test = zeros(params.number_drugs);
total_num_treated = zeros(params.number_drugs);
num_lost_resistance_to_test = zeros(params.number_drugs);

%for each antibiotic treated 
for drug = 1:params.number_drugs 
   total_prch(drug) = nnz(UTI_cases.PCR_sameday(:,drug)); 
   dates_index_drug =  find(UTI_cases.SamplingDate >= dates_to_use_start(drug) & UTI_cases.SamplingDate <= dates_to_use_end(drug));    
   
   %for each antibiotic suceptiblitly meased
   for drug_to_test = 1:params.number_drugs 
   dates_index_test =  find(UTI_cases.SamplingDate >= dates_to_use_start(drug_to_test) & UTI_cases.SamplingDate <= dates_to_use_end(drug_to_test));
   dates_index =  intersect(dates_index_test , dates_index_drug); 
   all_sensitive_test = (UTI_cases.SMP_Res(dates_index,drug_to_test) == 1 | UTI_cases.SMP_Res(dates_index,drug_to_test) == 2) & UTI_cases.hasdiag(dates_index) ;
   all_resistant_test = UTI_cases.SMP_Res(dates_index,drug_to_test) == 3  & UTI_cases.hasdiag(dates_index) ;
   all_currentres = ismember(UTI_cases.SMP_Res(dates_index,drug),[1 2]) & ismember(UTI_cases.SMP_Res(dates_index,drug_to_test),[1 2 3]);
   all_sensitive_nexttestres =  UTI_cases.next_res(dates_index,drug_to_test) == 1 | UTI_cases.SMP_Res(dates_index,drug_to_test) == 2 ;
   all_resistant_nexttestres =  UTI_cases.next_res(dates_index,drug_to_test) == 3 ;
   all_nextres = ismember(UTI_cases.next_res(dates_index,drug),[1 2 3]) & ismember(UTI_cases.next_res(dates_index,drug_to_test),[1 2 3]);
   gained_resistance_to_test = all_sensitive_test & all_resistant_nexttestres & UTI_cases.PCR_sameday(dates_index,drug);
   lost_resistance_to_test = all_resistant_test & all_sensitive_nexttestres & UTI_cases.PCR_sameday(dates_index,drug);
   total_num_treated(drug,drug_to_test) = nnz(UTI_cases.PCR_sameday(dates_index,drug) & all_nextres & all_currentres);
   num_gained_resistance_to_test(drug,drug_to_test) = nnz(gained_resistance_to_test & UTI_cases.treatfailure(dates_index) & all_nextres & all_currentres);
   num_lost_resistance_to_test(drug,drug_to_test) = nnz(lost_resistance_to_test & UTI_cases.treatfailure(dates_index) & all_nextres & all_currentres);
   end
end 

%% net change in resistance
rate_gained_resistance_to_test = num_gained_resistance_to_test./total_num_treated*100;
rate_lost_resistance_to_test = num_lost_resistance_to_test./total_num_treated*100;
rate_gained_resistance_to_test(rate_gained_resistance_to_test == Inf) = 0;
rate_lost_resistance_to_test(rate_lost_resistance_to_test == Inf) = 0;
net_change_res = rate_gained_resistance_to_test-rate_lost_resistance_to_test;
net_change_res(isnan(rate_lost_resistance_to_test)) = 0;
net_change_res = net_change_res(params.new_order, params.new_order);

n = 10;
lost_res_num = 15; % for colormap
gained_res_num = 55; % for colormap
total_num = gained_res_num + lost_res_num;
%make colormap 
cyan = [ linspace(0, 1, gained_res_num)' ones(gained_res_num,1) ones(gained_res_num,1)];
magenta = [  ones(gained_res_num,1) linspace(0, 1, gained_res_num)' ones(gained_res_num,1)];
magenta = flipud(magenta);
cyan_white_magenta = [cyan; magenta(1:lost_res_num,:)];
colormap(flipud(cyan_white_magenta))



%% chance of failure untreated
num_gained_resistance_to_test = zeros(params.number_drugs, 1);
total_num_treated = zeros(params.number_drugs,1);
num_lost_resistance_to_test = zeros(params.number_drugs,1);

for drug = 1:params.number_drugs 
   total_prch(drug) = nnz(UTI_cases.PCR_sameday(:,10));
   dates_index =  find(UTI_cases.SamplingDate >= dates_to_use_start(drug) & UTI_cases.SamplingDate <= dates_to_use_end(drug));
   all_sensitive_test = (UTI_cases.SMP_Res(dates_index,drug) == 1 | UTI_cases.SMP_Res(dates_index,drug) == 2) & UTI_cases.hasdiag(dates_index) ;
   all_resistant_test = UTI_cases.SMP_Res(dates_index,drug) == 3  & UTI_cases.hasdiag(dates_index) ;
   all_currentres = ismember(UTI_cases.SMP_Res(dates_index,drug),[1 2 3]);
   all_sensitive_nexttestres =  UTI_cases.next_res(dates_index,drug) == 1 | UTI_cases.SMP_Res(dates_index,drug) == 2 ;
   all_resistant_nexttestres =  UTI_cases.next_res(dates_index,drug) == 3 ;
   all_nextres = ismember(UTI_cases.next_res(dates_index,drug),[1 2 3]) & ismember(UTI_cases.next_res(dates_index,drug),[1 2 3]);
   gained_resistance_to_test = all_sensitive_test & all_resistant_nexttestres & UTI_cases.PCR_sameday(dates_index,10);
   lost_resistance_to_test = all_resistant_test & all_sensitive_nexttestres & UTI_cases.PCR_sameday(dates_index,10);
   total_num_treated(drug) = nnz(UTI_cases.PCR_sameday(dates_index,10) & all_nextres & all_currentres);
   num_gained_resistance_to_test(drug) = nnz(gained_resistance_to_test & UTI_cases.treatfailure(dates_index) & all_nextres & all_currentres);
   num_lost_resistance_to_test(drug) = nnz(lost_resistance_to_test & UTI_cases.treatfailure(dates_index) & all_nextres & all_currentres);

end 
%% net change in resistance
rate_gained_resistance_to_test = num_gained_resistance_to_test./total_num_treated*100;
rate_lost_resistance_to_test = num_lost_resistance_to_test./total_num_treated*100;
rate_gained_resistance_to_test(rate_gained_resistance_to_test == Inf) = 0;
rate_lost_resistance_to_test(rate_lost_resistance_to_test == Inf) = 0;
D3 = rate_gained_resistance_to_test-rate_lost_resistance_to_test;

%% make the plot

image([rot90(net_change_res) flipud(D3)]*n+lost_res_num )
xlabel('treatment antibiotic')
ylabel('susceptibiltiy of reccurent UTI')
colormap(flipud(cyan_white_magenta))
xticks([1:1:params.number_drugs+1])
set(gca,'xticklabel',[UTI_cases.SMP_Res_drug_names(params.new_order); 'Untreated'])
xtickangle(45);
yticks([1:1:params.number_drugs])
set(gca,'yticklabel',flipud(UTI_cases.SMP_Res_drug_names(params.new_order)))
cbh = colorbar ; %Create Colorbar
x=5;
cbh.Ticks = linspace(0, length(colormap), 8)+x ; 
ticknums = (linspace(-lost_res_num, gained_res_num, 8))/n+x/n;
cbh.TickLabels = num2cell(ticknums) ; 
set(get(cbh,'label'),'string',{'net % of treated cases which changed';...
    'resistance to focal drug'});
axis image

end

function [] = Patient_demographics_table_wounds2(wound_cases, params)
% funtion takes wound_cases struc and makes a table of patient demographics numbers and % 
% for all antibiotic treated cases. Writes table to the Tables directory. Table
% lists the following demographics:
names ={'Female'; 'Male';'Pregnant';'Age 0-9';'Age 10-19';...
    'Age 20-29';'Age 30-39';'Age 40-49';'Age 50-59';'Age 60-69';'Age 70-79';'Age 80-89';'Age >90';...
    'Prev SMPs 0-1'; 'Prev SMPs 2-4';'Prev SMPs 5+';};

%% only susceptibiltiy measured

for drug= 1:params.number_drugs
 relavant_suceptib_measured(:,drug) = wound_cases.PCR_sameday(:,drug) == 1  & wound_cases.SMP_Res(:, drug) ~= 0;
 relavant_suceptib_notmeasured(:,drug) = wound_cases.PCR_sameday(:,drug) == 1  & wound_cases.SMP_Res(:, drug) == 0;
 suceptib_matched(:,drug) =  wound_cases.PCR_sameday(:,drug) == 1  & ismember(wound_cases.SMP_Res(:, drug), [1 2]); 
 suceptib_mismatched(:,drug) =  wound_cases.PCR_sameday(:,drug) == 1  & ismember(wound_cases.SMP_Res(:, drug), 3); 
 treat_fail(:,drug) =  wound_cases.PCR_sameday(:,drug) == 1  & wound_cases.treatfailure; 
end

relavant_suceptib_measured_any = (sum(relavant_suceptib_measured,2));

for drug = 1:params.number_drugs
val1(drug) = nnz(relavant_suceptib_notmeasured(:,drug));
end

wound_cases.hasdiag = relavant_suceptib_measured_any; % just the susceptibiltiy measured cases
Demog_treated =  wound_cases.Demog(wound_cases.hasdiag & wound_cases.PCR_sameday_any ,:);
any_SRmeasurement_treated =  wound_cases.any_SRmeasurement( wound_cases.PCR_sameday_any,:,1); %+ ...

%%

predictive_numbers(1) = nnz(Demog_treated.Gender);
predictive_numbers(2) = nnz(~Demog_treated.Gender);
predictive_numbers(3) = nnz(Demog_treated.Preg);


for ii = 1:10
predictive_numbers(ii+3) = nnz(Demog_treated.Age(:,ii));    
end

prev_smp_range = [0 1 4 ];
any_SRmeasurement = max(any_SRmeasurement_treated,[],2);
predictive_numbers(14) = nnz(ismember(any_SRmeasurement	,prev_smp_range(1):prev_smp_range(2)));
predictive_numbers(15) = nnz(ismember(any_SRmeasurement	,prev_smp_range(2)+1:prev_smp_range(3)));
predictive_numbers(16) = nnz(any_SRmeasurement	> prev_smp_range(3));

predictive_per(1) = nnz(Demog_treated.Gender)/length(Demog_treated.Gender);
predictive_per(2) = nnz(~Demog_treated.Gender)/length(Demog_treated.Gender);
predictive_per(3) = nnz(Demog_treated.Preg)/length(Demog_treated.Gender);

for ii = 1:10
predictive_per(ii+3) = nnz(Demog_treated.Age(:,ii))/length(Demog_treated.Age(:,1)); 
end

predictive_per(14) = nnz(ismember(any_SRmeasurement	,prev_smp_range(1):prev_smp_range(2)))/length(any_SRmeasurement);
predictive_per(15) = nnz(ismember(any_SRmeasurement	,prev_smp_range(2)+1:prev_smp_range(3)))/length(any_SRmeasurement);
predictive_per(16) = nnz(any_SRmeasurement	> prev_smp_range(3))/length(any_SRmeasurement);
%% treatment failures
Demog_treated =  wound_cases.Demog(wound_cases.hasdiag & wound_cases.PCR_sameday_any & wound_cases.treatfailure,:);
any_SRmeasurement_treated =  wound_cases.any_SRmeasurement(wound_cases.hasdiag & wound_cases.PCR_sameday_any & wound_cases.treatfailure,:);

predictive_numbers_fails(1) = nnz(Demog_treated.Gender);
predictive_numbers_fails(2) = nnz(~Demog_treated.Gender);
predictive_numbers_fails(3) = nnz(Demog_treated.Preg);

for ii = 1:10 % age groups
predictive_numbers_fails(ii+3) = nnz(Demog_treated.Age(:,ii));    
end

any_SRmeasurement = max(any_SRmeasurement_treated,[],2);
predictive_numbers_fails(14) = nnz(ismember(any_SRmeasurement	,prev_smp_range(1):prev_smp_range(2)));
predictive_numbers_fails(15) = nnz(ismember(any_SRmeasurement	,prev_smp_range(2)+1:prev_smp_range(3)));
predictive_numbers_fails(16) = nnz(any_SRmeasurement > prev_smp_range(3));
predictive_per_fails(1) = nnz(Demog_treated.Gender)/length(Demog_treated.Gender);
predictive_per_fails(2) = nnz(~Demog_treated.Gender)/length(Demog_treated.Gender);
predictive_per_fails(3) = nnz(Demog_treated.Preg)/length(Demog_treated.Gender);

for ii = 1:10
predictive_per_fails(ii+3) = nnz(Demog_treated.Age(:,ii))/length(Demog_treated.Age(:,1)); 
end

predictive_per_fails(14) = nnz(ismember(any_SRmeasurement	,prev_smp_range(1):prev_smp_range(2)))/length(any_SRmeasurement);
predictive_per_fails(15) = nnz(ismember(any_SRmeasurement	,prev_smp_range(2)+1:prev_smp_range(3)))/length(any_SRmeasurement);
predictive_per_fails(16) = nnz(any_SRmeasurement	> prev_smp_range(3))/length(any_SRmeasurement);

%% make table

 Demographics_table_cases = table(names,predictive_numbers',round(predictive_per*100)'...
    ,predictive_numbers_fails',round(predictive_per_fails*100)');

filename = 'Tables/Patient_demographics_table_wounds.xlsx';
writetable(Demographics_table_cases,filename);

% sum(Demographics_table_cases.Var2(1:2))
% sum(Demographics_table_cases.Var2(5:14))
end
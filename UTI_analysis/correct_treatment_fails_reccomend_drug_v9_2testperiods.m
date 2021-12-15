function [] = correct_treatment_fails_reccomend_drug_v9_2testperiods(UTI_cases,params)
% function takes the UTI_case structure and optional params and treains and tests 
% regression to predict risk of gain of resistance for all susceptiblty-matched treated cases cases 
% model is trained in temporally seperated period, the cost adjusted in
% another, then tested in a third

figure
set(gcf,'color','w')
clf; set(gcf,'name','Fig. S14C seperate constrain vs test periods');
set(gcf,'units','centimeters','Position',[2 2 25 10]);

%demogs included
X_age = UTI_cases.Demog.Age;
X_age(:,6) = [];% reference age
zeroto39 = sum(X_age(:,1:4),2); % larger bin for younger ages since less data
eightyplus = sum(X_age(:,end-1:end),2);
X_demog = [zeroto39 X_age(:,5:end-2) eightyplus, UTI_cases.Demog.Gender, UTI_cases.Demog.Preg];
number_of_drugs = 7; % not including orflox as not measured in test period

sen_res = params.resistant_group; % gained resistance 

%% seperate training data and test data with seperate cost_ajust and test period
test2_date_start = datenum('2017-08-01') ;
test2_date_end = datenum('2018-05-01') ;

test_date_start = datenum('2018-05-01') ;
test_date_end = datenum('2019-07-01');

test_date_range = test_date_start:test_date_end;
test_data = ismember(UTI_cases.SamplingDate, test_date_range);
test_data2 = ismember(UTI_cases.SamplingDate, test2_date_start:test2_date_end);
fprintf('Testing using %.2f percent of the data\n' ,nnz(test_data)/length(test_data)*100)

%% Regression for StoR failure given any previous resistance to drug
clear coef

num_of_previous_R = UTI_cases.num_of_previous_R;
num_of_previous_S = UTI_cases.num_of_previous_S;

for drug = 1:number_of_drugs
SMP_to_use =  find((UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2)  & ~test_data & ~test_data2 & UTI_cases.hasdiag);
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), num_previous_R(drug_prch_temp==1), num_previous_S(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
[glm.B,glm.dev,glm.stats] = glmfit(X,Y,'binomial','logit') ;
coef(:,drug) = glm.B;
%nums(drug) = length(X);

end


%% test data:
clear P
P = nan(height(UTI_cases.Demog),number_of_drugs);
treatfails_gainedres_toprch = zeros(height(UTI_cases.Demog),number_of_drugs);
%patient was correctly prescribed a drug from the used drugs
was_correctly_prescribed = sum(UTI_cases.PCR_sameday(:,1:number_of_drugs) & ismember(UTI_cases.SMP_Res(:,1:number_of_drugs), [1 2]),2) > 0;
had_prev_resistance = nan(height(UTI_cases.Demog),number_of_drugs);

for drug = 1:number_of_drugs

SMP_to_use =  find(was_correctly_prescribed  & test_data2 & UTI_cases.hasdiag);

had_prev_resistance(was_correctly_prescribed  & test_data2 & UTI_cases.hasdiag & num_of_previous_R(:,drug) >0  ) = 1;   
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use, num_previous_R num_previous_S];
P(SMP_to_use, drug) = glmval(coef(:,drug),X,'logit') ;

%%
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
treatfails = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
treatfails_gainedres_toprch(SMP_to_use, drug) = treatfails & drug_prch_temp;

%% reassign any mismatched (or unknown Res) probabilities to 1
index = find(ismember(UTI_cases.SMP_Res(SMP_to_use,drug), [0 3]));
P(SMP_to_use(index), drug)  = 1.1;

end

treatfails_gainedres_toprch = logical(sum(treatfails_gainedres_toprch,2));

%% Test that reccomended which are incorectly prescribed are set to 1
fprintf('There minimum P for mismatched reccomendations is  %.1f \n' ,min(P(UTI_cases.SMP_Res(:,1:number_of_drugs) == 3)))


%% pick lowest P
clear M best_drug

withPdata_index = find(sum(isnan(P),2)== 0);
P_nonans = P(withPdata_index,:);
real_prch = UTI_cases.PCR_sameday(withPdata_index, 1:number_of_drugs);
real_treatfails_gainedres_toprch = treatfails_gainedres_toprch(withPdata_index);
real_res = UTI_cases.SMP_Res(withPdata_index, 1:number_of_drugs);

[M,best_drug_temp] = min(P_nonans,[],2);

best_drug = zeros(size(P_nonans,1),number_of_drugs);
for ii = 1:number_of_drugs
best_drug(best_drug_temp == ii,ii) = 1;
end


%% Test that data all have a prch amoung the chosen drugs 
fprintf('There are %.0f samples without a PRCH (should be zero)\n' ,(length(real_prch) - nnz(sum(real_prch,2)>0)))

% Test that there are no actually inccorectly prescribed samples 
incorrectly_presribed = nnz(sum(ismember(real_res, 3) & real_prch,2));
fprintf('There are %.0f samples incorrectly prescribed (should be zero)\n' ,incorrectly_presribed)

correctly_presribed = nnz(sum(ismember(real_res, [1 2]) & real_prch,2));
fprintf('%.0f out of %.0f samples are actuall correctly prescribed (should be all)\n' ,correctly_presribed, length(real_res))

reccomended_correctly_presribed = logical(sum(ismember(real_res, [1 2]) & best_drug,2));
fprintf('%.0f out of %.0f BEST RECCOMENDED samples are correctly prescribed (should be all)\n' ,nnz(reccomended_correctly_presribed), length(reccomended_correctly_presribed))


%% are the NON COST ADJUSTED prescribed and best reccomended the same?
prescribed_pred_matched = logical(sum(best_drug & real_prch,2));

fprintf('Mean rate of S->R fail for pysician/recommended (no cost adj) MATCHED is: %.2f %%\n' ,mean(real_treatfails_gainedres_toprch(prescribed_pred_matched))*100)
fprintf('Mean rate of S->R fail for pysician/recommended (no cost adj) NOT MATCHED is: %.2f %%\n' ,mean(real_treatfails_gainedres_toprch(~prescribed_pred_matched))*100)

mean(real_treatfails_gainedres_toprch)
%mean(real_fails2(~prescribed_pred_matched))

%% was there alternative available
num_w_alternative = nnz(best_drug & ~real_prch)/ (nnz(best_drug & real_prch) + nnz(best_drug & ~real_prch))*100;
disp(['in ' num2str(num_w_alternative) ' of cases an alternative (non cost adj) drug w/ lower P was avavilable']); 
%% cost adjusted to give equal numbers of reccomended drugs to actual purchased drugs 
rng('default'); rng(1);
P_nonans = P_nonans+rand(size(P_nonans))/100000; % to help with the cost adjusting

Fr0= zeros(1, number_of_drugs);
for ii = 1:number_of_drugs
Fr0(ii) = nnz(real_prch(:,ii)); 
end
Fr0 = Fr0./size(real_prch,1);

iter_params = [1 0.9 0.9 1e4];
 % params = [0.1, 0.9, 0.95, 1e4] ;
d = size(P_nonans,2) ;
n = size(P_nonans,1) ;

tol = iter_params(2) / n ;
relax = iter_params(3) ;
niter = iter_params(4) ;

dcost = iter_params(1) * ones(1,d) ;
cost = zeros(1,d) ;
exitflag = 1 ;
Fr = Fr0 ;

for i = 1:niter
    cost = cost + dcost.*sign(Fr-Fr0) ;
    cost = cost - mean(cost) ;
    dcost = dcost*relax ;
    
    [p,treat] = min(P_nonans+cost,[],2) ;
    treat(p>1) = 0 ;
    Fr = histc(treat,1:d)' ./ n ;
    err = max(abs(Fr-Fr0)) ;
    %Fr
    if err<tol; exitflag = 0 ; break; end
end


if err
    fprintf('error in number of samples=\n')
    disp((Fr-Fr0)*n) ;
end

%% test data:
clear P
P = nan(height(UTI_cases.Demog),number_of_drugs);
treatfails_gainedres_toprch = zeros(height(UTI_cases.Demog),number_of_drugs);
%patient was correctly prescribed a drug from the used drugs
was_correctly_prescribed = sum(UTI_cases.PCR_sameday(:,1:number_of_drugs) & ismember(UTI_cases.SMP_Res(:,1:number_of_drugs), [1 2]),2) > 0;

had_prev_resistance = nan(height(UTI_cases.Demog),number_of_drugs);

for drug = 1:number_of_drugs

SMP_to_use =  find(was_correctly_prescribed  & test_data & UTI_cases.hasdiag);

had_prev_resistance(was_correctly_prescribed  & test_data & UTI_cases.hasdiag & num_of_previous_R(:,drug) >0  ) = 1;   
any_previous_R = num_of_previous_R(SMP_to_use ,drug)>0 ;
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);

X_demog_to_use = X_demog(SMP_to_use,:);

X = [X_demog_to_use, num_previous_R num_previous_S];

P(SMP_to_use, drug) = glmval(coef(:,drug),X,'logit') ;

%%
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
treatfails = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
treatfails_gainedres_toprch(SMP_to_use, drug) = treatfails & drug_prch_temp;

%% reassign any mismatched (or unknown Res) probabilities to 1
index = find(ismember(UTI_cases.SMP_Res(SMP_to_use,drug), [0 3]));
P(SMP_to_use(index), drug)  = 1.1;

end

treatfails_gainedres_toprch = logical(sum(treatfails_gainedres_toprch,2));

%% Test that reccomended which are incorectly prescribed are set to 1
fprintf('There minimum P for mismatched reccomendations is  %.1f \n' ,min(P(UTI_cases.SMP_Res(:,1:number_of_drugs) == 3)))


%% pick lowest P
clear M best_drug

withPdata_index = find(sum(isnan(P),2)== 0);
P_nonans = P(withPdata_index,:);
real_prch = UTI_cases.PCR_sameday(withPdata_index, 1:number_of_drugs);
real_treatfails_gainedres_toprch = treatfails_gainedres_toprch(withPdata_index);
real_res = UTI_cases.SMP_Res(withPdata_index, 1:number_of_drugs);

[M,best_drug_temp] = min(P_nonans,[],2);

best_drug = zeros(size(P_nonans,1),number_of_drugs);
for ii = 1:number_of_drugs
best_drug(best_drug_temp == ii,ii) = 1;
end


%% Test that data all have a prch amoung the chosen drugs 
fprintf('There are %.0f samples without a PRCH (should be zero)\n' ,(length(real_prch) - nnz(sum(real_prch,2)>0)))

% Test that there are no actually inccorectly prescribed samples 
incorrectly_presribed = nnz(sum(ismember(real_res, 3) & real_prch,2));
fprintf('There are %.0f samples incorrectly prescribed (should be zero)\n' ,incorrectly_presribed)

correctly_presribed = nnz(sum(ismember(real_res, [1 2]) & real_prch,2));
fprintf('%.0f out of %.0f samples are actuall correctly prescribed (should be all)\n' ,correctly_presribed, length(real_res))

reccomended_correctly_presribed = logical(sum(ismember(real_res, [1 2]) & best_drug,2));
fprintf('%.0f out of %.0f BEST RECCOMENDED samples are correctly prescribed (should be all)\n' ,nnz(reccomended_correctly_presribed), length(reccomended_correctly_presribed))


%% are the NON COST ADJUSTED prescribed and best reccomended the same?
prescribed_pred_matched = logical(sum(best_drug & real_prch,2));

fprintf('Mean rate of S->R fail for pysician/recommended (no cost adj) MATCHED is: %.2f %%\n' ,mean(real_treatfails_gainedres_toprch(prescribed_pred_matched))*100)
fprintf('Mean rate of S->R fail for pysician/recommended (no cost adj) NOT MATCHED is: %.2f %%\n' ,mean(real_treatfails_gainedres_toprch(~prescribed_pred_matched))*100)

mean(real_treatfails_gainedres_toprch)
%mean(real_fails2(~prescribed_pred_matched))

%% was there alternative available
num_w_alternative = nnz(best_drug & ~real_prch)/ (nnz(best_drug & real_prch) + nnz(best_drug & ~real_prch))*100;
disp(['in ' num2str(num_w_alternative) ' of cases an alternative (non cost adj) drug w/ lower P was avavilable']); 
    
    %%
 %add previous cost to the new Ps   
[p,treat] = min(P_nonans+cost,[],2) ;

%check how good the prev cost on the new data is
Fr0= zeros(1, number_of_drugs);
for ii = 1:number_of_drugs
Fr0(ii) = nnz(real_prch(:,ii)); 
end
Fr0 = Fr0./size(real_prch,1);
d = size(P_nonans,2) ;
n = size(P_nonans,1) 
Fr = histc(treat,1:d)' ./ n ;
fprintf('error in number of samples=\n')
    disp((Fr-Fr0)*100) ;
    mean(abs((Fr-Fr0)*100)) 
best_drug_costadjust = zeros(size(treat,1),number_of_drugs);
for ii = 1:number_of_drugs
best_drug_costadjust(treat == ii,ii) = 1;

end

%random but correctly prescribed
random_drug = zeros(size(treat,1),number_of_drugs);
num = 1:7;
aa = 0;
for ii = 1:size(treat,1)
numstemp = num(P_nonans(ii,:) < 1);
if ~isempty(numstemp)
random_drug(ii,numstemp(randi(length(numstemp)))) = 1;
else
    numstemp = 1:7;
    random_drug(ii,numstemp(randi(length(numstemp)))) = 0;
    aa =aa+1;
end
end

had_prev_resistance_nonans = num_of_previous_R(withPdata_index,1:7);
incorrectly_prescribed = find(sum(had_prev_resistance_nonans & real_prch,2)>0);
% size(had_prev_resistance_nonans)
% size(real_prch)
random_drug_noprev = real_prch;
num = 1:7;
for ii = 1:length(incorrectly_prescribed)
    random_drug_noprev(incorrectly_prescribed(ii),1:7)= zeros(1,7);
numstemp = num(had_prev_resistance_nonans(incorrectly_prescribed(ii),1:7) == 0 & ismember(real_res(incorrectly_prescribed(ii),1:7), [1 2]));
if ~isempty(numstemp)
random_drug_noprev(incorrectly_prescribed(ii),numstemp(randi(length(numstemp)))) = 1;
end
end

%%
numprch = zeros(1,number_of_drugs);
numreccost = zeros(1,number_of_drugs);
numrec = zeros(1,number_of_drugs);
for ii = 1:number_of_drugs 
numprch(ii)= nnz(real_prch(:,ii));
numreccost(ii)= nnz(best_drug_costadjust(:,ii));
numrec(ii)= nnz(best_drug(:,ii));
end

%% cost adjusted reccomended only  correctly prescribed

real = mean(P_nonans(real_prch ==1));
best_cost = mean(P_nonans(best_drug_costadjust ==1));
best = mean(P_nonans(best_drug ==1));
random_noprev = mean(P_nonans(random_drug_noprev ==1));

realPs = P_nonans(real_prch ==1);
best_costPs = P_nonans(best_drug_costadjust ==1);
bestPs = P_nonans(best_drug ==1);
random_noprevPs = P_nonans(random_drug_noprev ==1);
num_bootstraps = 1000;

real_boot = bootstrp(num_bootstraps,@mean,realPs);
best_cost_boot = bootstrp(num_bootstraps,@mean,best_costPs);
best_boot = bootstrp(num_bootstraps,@mean,bestPs);
random_noprev_boot = bootstrp(num_bootstraps,@mean,random_noprevPs);

random_noprev_CI = prctile(random_noprev_boot,[2.5 97.5]) ;
real_CI = prctile(real_boot,[2.5 97.5]) ;
best_CI = prctile(best_boot,[2.5 97.5]) ;
best_cost_CI = prctile(best_cost_boot,[2.5 97.5]) ;

subplot(1,2,2)
b1= bar([mean(real_boot) mean(random_noprev_boot)  mean(best_cost_boot) mean(best_boot)]);
hold on
%errorbar(b1(1).XData,b1(1).YData,[std(random_boot) std(real_boot)  std(best_cost_boot) std(best_boot)],'k', 'LineStyle','none')
errorbar(b1(1).XData,b1(1).YData,[ mean(real_boot)-real_CI(1) mean(random_noprev_boot)-random_noprev_CI(1) mean(best_cost_boot)-best_cost_CI(1) mean(best_boot)-best_CI(1)  ],...
    [ -mean(real_boot)+real_CI(2) -mean(random_noprev_boot)+random_noprev_CI(1) -mean(best_cost_boot)+best_cost_CI(2) -mean(best_boot)+best_CI(2)  ],'k', 'LineStyle','none')
xticklabels({'Physician','Random no prevR','Recomended cost adj','Recommended'})
xtickangle( 45 )
ylabel('predicted probability of aquiring resistance')
plot(1:4,mean(real_treatfails_gainedres_toprch)*ones(1,4),'k--');

%% which drug best vs physician prescibed

number_prescribed_reccomended = zeros(3,size(real_prch,2));
number_prescribed_reccomended(1,:) = sum(real_prch);
number_prescribed_reccomended(2,:) = sum(best_drug_costadjust);
number_prescribed_reccomended(3,:) = sum(best_drug);


subplot(1,2,1)
bar(number_prescribed_reccomended')
legend({'physician prescribed drugs' ; 'ML, recommended cost adjusted' ; 'ML, recommended'})
xticklabels(UTI_cases.SMP_Res_drug_names)
xtickangle(45)
ylabel('number of prescriptions/ reccomendations')

end
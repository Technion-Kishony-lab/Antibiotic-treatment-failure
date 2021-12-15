function [] = reccomend_drugs(wound_cases, params)
% function takes the UTI_case structure and optional params and treains and tests 
% regressoin model to predict risk of gain of resistance for all susceptibility-matched treated cases cases 

X_age = wound_cases.Demog.Age; X_age(:,6) = []; % reference age
zeroto39 = sum(X_age(:,1:2),2); % larger bin for lower age groups with less data
eightyplus = sum(X_age(:,end-1:end),2);
X_demog = [zeroto39 X_age(:,3:end) , wound_cases.Demog.Gender ];

%% seperate training data and test data
test_data = ismember(wound_cases.SamplingDate, 736411:737382);
fprintf('Testing using %.2f percent of the data\n' ,nnz(test_data)/length(test_data)*100)

%% Regression for StoR failure given any previous resistance to drug
clear coef
num_of_previous_R = wound_cases.num_of_previous_R;
num_of_previous_S = wound_cases.num_of_previous_S;
sen_res = params.resistant_group; % gain of resistance
warning('off','stats:glmfit:IterationLimit')

for drug = 1:params.number_drugs
SMP_to_use =  find((wound_cases.SMP_Res(:,drug)== 1 | wound_cases.SMP_Res(:,drug)== 2)  & ~test_data ); 
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);
treatfailure = wound_cases.treatfailure(SMP_to_use) ==1 & ismember(wound_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = wound_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), num_previous_R(drug_prch_temp==1), num_previous_S(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
[glm.B,glm.dev,glm.stats] = glmfit(X,Y,'binomial','logit') ;
coef(:,drug) = glm.B;
%AUC
scores = glmval(glm.B,X,'logit') ;
[X,Y,T,AUC] = perfcurve(Y,scores,1);
AUC_all(drug) = AUC;
 subplot(4,2,drug)
 plot(X,Y);
 xlabel('False positive rate') ;
 ylabel('True positive rate');
 title([wound_cases.SMP_Res_drug_names{drug} '  AUC:' num2str(AUC)]);
end

%% test data
clear P
P = nan(height(wound_cases.Demog),params.number_drugs);
treatfails_gainedres_toprch = zeros(height(wound_cases.Demog),params.number_drugs);
to_use = false(height(wound_cases.Demog),params.number_drugs);
%patient was correctly prescribed a drug from the used drugs
was_correctly_prescribed = sum(wound_cases.PCR_sameday(:,1:params.number_drugs) & ismember(wound_cases.SMP_Res(:,1:params.number_drugs), [1 2]),2) > 0;

had_prev_resistance = nan(height(wound_cases.Demog),params.number_drugs);

for drug = 1:params.number_drugs
SMP_to_use =  find(was_correctly_prescribed  & test_data );
had_prev_resistance(was_correctly_prescribed  & test_data  & num_of_previous_R(:,drug) >0  ) = 1;   
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use, num_previous_R num_previous_S];
P(SMP_to_use, drug) = glmval(coef(:,drug),X,'logit') ;
%
drug_prch_temp = wound_cases.PCR_sameday(SMP_to_use ,drug);
treatfails = wound_cases.treatfailure(SMP_to_use) ==1 & ismember(wound_cases.next_res(SMP_to_use, drug),sen_res) ;
treatfails_gainedres_toprch(SMP_to_use, drug) = treatfails & drug_prch_temp;
num_gained(drug) = nnz(treatfails(drug_prch_temp));
num_treated(drug) = nnz(drug_prch_temp);
% reassign any mismatched (or unknown Res) probabilities to 1
index = find(ismember(wound_cases.SMP_Res(SMP_to_use,drug), [0 3]));
P(SMP_to_use(index), drug)  = 1.1; % mismatched treatments are not used 

X_test = [X_demog_to_use(drug_prch_temp==1,:), num_previous_R(drug_prch_temp==1), num_previous_S(drug_prch_temp==1)];
treatfailure = wound_cases.treatfailure(SMP_to_use) ==1 & ismember(wound_cases.next_res(SMP_to_use, drug),sen_res) ;
Y = treatfailure(drug_prch_temp==1);

scores = glmval(coef(:,drug),X_test,'logit') ;
[X_test,Y,T,AUC2] = perfcurve(Y,scores,1);
AUC_all2(drug) = AUC2;
 subplot(4,2,drug)
 hold on
 plot(X_test,Y,'r');
 xlabel('False positive rate') ;
 ylabel('True positive rate');
 title([wound_cases.SMP_Res_drug_names{drug} '  AUC train:' num2str(AUC_all(drug)) 'AUC test:' num2str(AUC_all2(drug))]);
end

sum(num_gained)/sum(num_treated)
treatfails_gainedres_toprch = logical(sum(treatfails_gainedres_toprch,2));

%% pick antibiotic with lowest P
clear M best_drug

withPdata_index = find(sum(isnan(P),2)== 0);
P_nonans = P(withPdata_index,:);
real_prch = wound_cases.PCR_sameday(withPdata_index, 1:params.number_drugs);
real_treatfails_gainedres_toprch = treatfails_gainedres_toprch(withPdata_index(1:end-100));
real_res = wound_cases.SMP_Res(withPdata_index, 1:params.number_drugs);

[M,best_drug_temp] = min(P_nonans,[],2);

best_drug = zeros(size(P_nonans,1),params.number_drugs);
for ii = 1:params.number_drugs
best_drug(best_drug_temp == ii,ii) = 1;
end

%% was there alternative available
num_w_alternative = nnz(best_drug & ~real_prch)/ (nnz(best_drug & real_prch) + nnz(best_drug & ~real_prch))*100;
disp(['in ' num2str(num_w_alternative) ' of cases an alternative (non cost adj) drug w/ lower P was avavilable']); 
%% cost adjusted to give equal numbers of reccomended drugs to actual purchased drugs 
rng('default'); rng(1);
P_nonans = P_nonans+rand(size(P_nonans))/100000; % to help with the cost adjusting
Fr0= zeros(1, params.number_drugs);
for ii = 1:params.number_drugs
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
%% pick antibiotic with lowest P after cost-adjusting P  
best_drug_costadjust = zeros(size(treat,1),params.number_drugs);
for ii = 1:params.number_drugs
best_drug_costadjust(treat == ii,ii) = 1;

end

%% avoid antibiotics to which patient had previous resistant infection
% keep actual treatment if no prev resistance infection
had_prev_resistance_nonans = num_of_previous_R(withPdata_index,1:5);
incorrectly_prescribed = find(sum(had_prev_resistance_nonans & real_prch,2)>0);
random_drug_noprev = real_prch;
num = 1:5;
for ii = 1:length(incorrectly_prescribed)
    random_drug_noprev(incorrectly_prescribed(ii),1:5)= zeros(1,5);
numstemp = num(had_prev_resistance_nonans(incorrectly_prescribed(ii),1:5) == 0 & ismember(real_res(incorrectly_prescribed(ii),1:5), [1 2]));
if ~isempty(numstemp)
random_drug_noprev(incorrectly_prescribed(ii),numstemp(randi(length(numstemp)))) = 1;
end
end

%% bootstrap data and plot predicted risk
reccomended_cost_correctly_presribed = logical(sum(ismember(real_res, [1 2]) & best_drug_costadjust,2));
fprintf('%.0f out of %.0f COST ADJUSTED BEST RECCOMENDED samples are correctly prescribed (should be all)\n' ,nnz(reccomended_cost_correctly_presribed), length(reccomended_cost_correctly_presribed))

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

figure;
set(gcf,'color','w', 'name','Fig. 4E', 'units','centimeters','Position',[1 1 10 10]);
b1= bar([mean(real_boot) mean(random_noprev_boot)  mean(best_cost_boot) mean(best_boot)]);
hold on
errorbar(b1(1).XData,b1(1).YData,...
    [ mean(real_boot)-real_CI(1) mean(random_noprev_boot)-random_noprev_CI(1) mean(best_cost_boot)-best_cost_CI(1) mean(best_boot)-best_CI(1)  ],...
    [ -mean(real_boot)+real_CI(2) -mean(random_noprev_boot)+random_noprev_CI(2) -mean(best_cost_boot)+best_cost_CI(2) -mean(best_boot)+best_CI(2)  ],'k', 'LineStyle','none')
xticklabels({'Physician','Random no prevR','Recomended cost adj','Recommended'})
xtickangle( 45 )
ylabel('predicted probability of aquiring resistance')
plot(1:4,mean(real_treatfails_gainedres_toprch)*ones(1,4),'k--');

%% number of each antibiotic prescribed in test period for 
% different prescription methods 
number_prescribed_reccomended = zeros(3,size(real_prch,2));
number_prescribed_reccomended(1,:) = sum(real_prch);
number_prescribed_reccomended(2,:) = sum(best_drug_costadjust);
number_prescribed_reccomended(3,:) = sum(best_drug);

figure
set(gcf,'color','w')
clf; set(gcf,'name','Fig. S14 Number of prescribed and ML recommended antibiotics');
bar(number_prescribed_reccomended')
legend({'physician prescribed drugs' ; 'ML, recommended cost adjusted' ; 'ML, recommended'})
xticklabels(wound_cases.SMP_Res_drug_names)
xtickangle(45)
ylabel('number of prescriptions/ reccomendations')

%% errors between models
dprb = real_boot-best_boot ;
sd = std(dprb) ; mn = mean(dprb) ; pval = 1-normcdf(mn/sd) 
dprb = real_boot-best_cost_boot ;
sd = std(dprb) ; mn = mean(dprb) ; pval = 1-normcdf(mn/sd) 
dprb = real_boot-random_noprev_boot ;
sd = std(dprb) ; mn = mean(dprb) ; pval = 1-normcdf(mn/sd)

%% Do regression for any failure (S-S & S-R) and look at the predicted risk for the 
% antibiotic recomendations trained to miniize StoR failures
clear coef
sen_res = [1 2 3];
for drug = 1:params.number_drugs
SMP_to_use =  find((wound_cases.SMP_Res(:,drug)== 1 | wound_cases.SMP_Res(:,drug)== 2)  & ~test_data );
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);
treatfailure = wound_cases.treatfailure(SMP_to_use) ==1 & ismember(wound_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = wound_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), num_previous_R(drug_prch_temp==1), num_previous_S(drug_prch_temp==1)];
%X = [X_demog_to_use(drug_prch_temp==1,:), any_previous_R(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
[glm.B,glm.dev,glm.stats] = glmfit(X,Y,'binomial','logit') ;
coef(:,drug) = glm.B;
end

%% test data:
clear P
P = nan(height(wound_cases.Demog),params.number_drugs);
treatfails_gainedres_toprch = zeros(height(wound_cases.Demog),params.number_drugs);
treatfails_all = zeros(height(wound_cases.Demog),params.number_drugs);
%patient was correctly prescribed a drug from the used drugs
was_correctly_prescribed = sum(wound_cases.PCR_sameday(:,1:params.number_drugs) & ismember(wound_cases.SMP_Res(:,1:params.number_drugs), [1 2]),2) > 0;

for drug = 1:params.number_drugs

SMP_to_use =  find(was_correctly_prescribed  & test_data ); 
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);
% treatfailure = SMP.treatfailure(SMP_to_use) ==1 & ismember(next_res(SMP_to_use, drug),sen_res) ;
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use, num_previous_R num_previous_S];
P(SMP_to_use, drug) = glmval(coef(:,drug),X,'logit') ;
drug_prch_temp = wound_cases.PCR_sameday(SMP_to_use ,drug);
treatfails = wound_cases.treatfailure(SMP_to_use) ==1 & ismember(wound_cases.next_res(SMP_to_use, drug),sen_res) ;

% reassign any mismatched (or unknown Res) probabilities to 1
%index = find(SMP_Res(SMP_to_use,drug)== 3);
index = find(ismember(wound_cases.SMP_Res(SMP_to_use,drug), [0 3]));
P(SMP_to_use(index), drug)  = 1.1;

end

withPdata_index_SS = find(sum(isnan(P),2)== 0);
P_nonans_SS = P(withPdata_index_SS,:);
best_cost_SS_boot = bootstrp(num_bootstraps,@mean,P_nonans_SS(best_drug_costadjust ==1));
best_SS_boot = bootstrp(num_bootstraps,@mean,P_nonans_SS(best_drug ==1));
real_SS_boot = bootstrp(num_bootstraps,@mean,P_nonans_SS(real_prch ==1));
%
figure;
set(gcf,'color','w', 'name','Fig. S17', 'units','centimeters','Position',[1 1 10 10]);
b1= bar([mean(real_SS_boot) mean(best_SS_boot) mean(best_cost_SS_boot) ]);
hold on
errorbar(b1(1).XData,b1(1).YData,[std(real_SS_boot) std(best_SS_boot) std(best_cost_SS_boot)  ],'k', 'LineStyle','none')
xticklabels({'All rec real','All rec best','All rec cost'})
xtickangle( 45 )
ylabel('predicted probability of aquiring resistance')
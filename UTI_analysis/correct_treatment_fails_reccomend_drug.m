function [] = correct_treatment_fails_reccomend_drug(UTI_cases,params)
% function takes the UTI_case structure and optional params and treains and tests 
% regressoin model to predict risk of gain of resistance for all susceptibility-matched treated cases cases 


X_age = UTI_cases.Demog.Age;
X_age(:,6) = []; % reference age
X_demog = [X_age, UTI_cases.Demog.Gender, UTI_cases.Demog.Preg UTI_cases.Demog.any_prev_cath];
number_of_drugs = 7; % not including orflox as not measured in test period
sen_res = params.resistant_group;
%% seperate training data and test data
test_date_start = datenum('2018-05-01') ;
test_date_end = datenum('2019-07-01') ;
test_date_range = test_date_start:test_date_end;
test_data = ismember(UTI_cases.SamplingDate, test_date_range);
fprintf('Testing using %.2f percent of the data\n' ,nnz(test_data)/length(test_data)*100)

%% Regression for StoR failure given any previous resistance to drug
clear coef
warning('off','stats:glmfit:IterationLimit');
any_SRmeasurement = UTI_cases.any_SRmeasurement;
num_of_previous_R = UTI_cases.num_of_previous_R;
num_of_previous_S = UTI_cases.num_of_previous_S;

figure;
set(gcf,'color','w', 'name','Fig. S12', 'units','centimeters','Position',[1 1 15 15]);

%loop over all antibiotics 
for drug = 1:number_of_drugs
%use cases sentivie to antibiotic     
SMP_to_use =  find( ismember(UTI_cases.SMP_Res(:,drug), params.sensitive_group)  & ~test_data & UTI_cases.hasdiag);
%number past infections sen or res to the antibiotic 
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);
%gained resistance treatment failure
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
%treated with antibiotic
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), num_previous_R(drug_prch_temp==1), num_previous_S(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
[glm.B,glm.dev,glm.stats] = glmfit(X,Y,'binomial','logit') ;
coef(:,drug) = glm.B;
% plots the AUC for training period
scores = glmval(glm.B,X,'logit') ;
[X,Y,T,AUC] = perfcurve(Y,scores,1);
AUC_all(drug) = AUC;
subplot(4,2,drug)
plot(X,Y);
xlabel('False positive rate') ;
ylabel('True positive rate');
title([UTI_cases.SMP_Res_drug_names{drug} '  AUC:' num2str(AUC)]);

end

%% test data:
clear P
P = nan(height(UTI_cases.Demog),number_of_drugs);
treatfails_gainedres_toprch = zeros(height(UTI_cases.Demog),number_of_drugs);
%patient was correctly prescribed a drug from the used drugs
was_correctly_prescribed = sum(UTI_cases.PCR_sameday(:,1:number_of_drugs) & ismember(UTI_cases.SMP_Res(:,1:number_of_drugs), params.sensitive_group),2) > 0;

for drug = 1:number_of_drugs

SMP_to_use =  find(was_correctly_prescribed & test_data & UTI_cases.hasdiag);
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use, num_previous_R num_previous_S];
P(SMP_to_use, drug) = glmval(coef(:,drug),X,'logit') ;
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
treatfails = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
treatfails_gainedres_toprch(SMP_to_use, drug) = treatfails & drug_prch_temp;
% reassign any mismatched (or unknown Res) probabilities to > 1
index = find(ismember(UTI_cases.SMP_Res(SMP_to_use,drug), [0 3]));
P(SMP_to_use(index), drug)  = 1.1; % mismatched treatments are not used 
X_test = [X_demog_to_use(drug_prch_temp==1,:), num_previous_R(drug_prch_temp==1), num_previous_S(drug_prch_temp==1)];
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
Y = treatfailure(drug_prch_temp==1);
scores = glmval(coef(:,drug),X_test,'logit') ;
% plots the AUC for test period
[X_test,Y,T,AUC2] = perfcurve(Y,scores,1);
AUC_all2(drug) = AUC2;
subplot(4,2,drug)
hold on
plot(X_test,Y,'r');
xlabel('False positive rate') ;
ylabel('True positive rate');
title([UTI_cases.SMP_Res_drug_names{drug} '  AUC train:' num2str(AUC_all(drug)) 'AUC test:' num2str(AUC_all2(drug))]);
end
treatfails_gainedres_toprch = logical(sum(treatfails_gainedres_toprch,2));

%% pick antibiotic with lowest P
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
incorrectly_presribed = nnz(sum(ismember(real_res, params.resistant_group) & real_prch,2));
fprintf('There are %.0f samples incorrectly prescribed (should be zero)\n' ,incorrectly_presribed)

correctly_presribed = nnz(sum(ismember(real_res, params.sensitive_group) & real_prch,2));
fprintf('%.0f out of %.0f samples are actuall correctly prescribed (should be all)\n' ,correctly_presribed, length(real_res))

reccomended_correctly_presribed = logical(sum(ismember(real_res, params.sensitive_group) & best_drug,2));
fprintf('%.0f out of %.0f BEST RECCOMENDED samples are correctly prescribed (should be all)\n' ,nnz(reccomended_correctly_presribed), length(reccomended_correctly_presribed))

%% are the NON COST ADJUSTED prescribed and best reccomended the same?
prescribed_pred_matched = logical(sum(best_drug & real_prch,2));
fprintf('Mean rate of S->R fail for pysician/recommended (no cost adj) MATCHED is: %.2f %%\n' ,mean(real_treatfails_gainedres_toprch(prescribed_pred_matched))*100)
fprintf('Mean rate of S->R fail for pysician/recommended (no cost adj) NOT MATCHED is: %.2f %%\n' ,mean(real_treatfails_gainedres_toprch(~prescribed_pred_matched))*100)

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
d = size(P_nonans,2) ;
n = size(P_nonans,1) ;
tol = iter_params(2) / n ;
relax = iter_params(3) ;
niter = iter_params(4) ;
dcost = iter_params(1) * ones(1,d) ;
cost = zeros(1,d) ;
exitflag = 1 ;
Fr = Fr0 ;
%iterative cost ajust P
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
best_drug_costadjust = zeros(size(treat,1),number_of_drugs);
for ii = 1:number_of_drugs
best_drug_costadjust(treat == ii,ii) = 1;
end

%% avoid antibiotics to which patient had previous resistant infection
% keep actual treatment if no prev resistance infection
had_prev_resistance_nonans = num_of_previous_R(withPdata_index,1:number_of_drugs);
incorrectly_prescribed = find(sum(had_prev_resistance_nonans & real_prch,2)>0);
random_drug_noprev = real_prch;
num = 1:number_of_drugs;
for ii = 1:length(incorrectly_prescribed)
    random_drug_noprev(incorrectly_prescribed(ii),1:number_of_drugs)= zeros(1,number_of_drugs);
numstemp = num(had_prev_resistance_nonans(incorrectly_prescribed(ii),1:number_of_drugs) == 0 & ismember(real_res(incorrectly_prescribed(ii),1:number_of_drugs), params.sensitive_group));
if ~isempty(numstemp)
random_drug_noprev(incorrectly_prescribed(ii),numstemp(randi(length(numstemp)))) = 1;
end
end


%% bootstrap data and plot predicted risk
reccomended_cost_correctly_presribed = logical(sum(ismember(real_res, params.sensitive_group) & best_drug_costadjust,2));
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
%errorbar(b1(1).XData,b1(1).YData,[std(random_boot) std(real_boot)  std(best_cost_boot) std(best_boot)],'k', 'LineStyle','none')
errorbar(b1(1).XData,b1(1).YData,[ mean(real_boot)-real_CI(1) mean(random_noprev_boot)-random_noprev_CI(1) mean(best_cost_boot)-best_cost_CI(1) mean(best_boot)-best_CI(1)  ],...
    [ -mean(real_boot)+real_CI(2) -mean(random_noprev_boot)+random_noprev_CI(1) -mean(best_cost_boot)+best_cost_CI(2) -mean(best_boot)+best_CI(2)  ],'k', 'LineStyle','none')
xticklabels({'Physician','Random no prevR','Recomended cost adj','Recommended'})
xtickangle( 45 )
ylabel('predicted probability of aquiring resistance')
plot(1:4,mean(real_treatfails_gainedres_toprch)*ones(1,4),'k--');

%% errors between models
dprb = real_boot-best_boot ;
sd = std(dprb) ; mn = mean(dprb) ; pval = 1-normcdf(mn/sd) 
dprb = real_boot-best_cost_boot ;
sd = std(dprb) ; mn = mean(dprb) ; pval = 1-normcdf(mn/sd) 

%% plot bar chart of gained-resistance rate for cases split by if the treatment reccomended or not

clear x
below_thresh = zeros(size(P_nonans));
for drug = 1:number_of_drugs
threshold(drug) = prctile(P_nonans(P_nonans(:,drug)<1,drug),85); % binarize into top 15% bottom 85%. See paper
below_thresh(:,drug) = P_nonans(:,drug) < threshold(drug);
below_thresh(P_nonans(:,drug) > 1,drug) = NaN;
end

above_thresh = below_thresh == 0;
below_thresh = below_thresh == 1;

was_reccomended = logical(real_prch) & ismember(real_res, [ 1 2]) & below_thresh;
was_not_reccomended = logical(real_prch) & ismember(real_res, [ 1 2]) & above_thresh;

p = zeros(1,number_of_drugs);
for drug = 1:number_of_drugs
   x(1,1,drug) =  nnz(real_treatfails_gainedres_toprch(was_reccomended(:,drug)));
   x(2,1,drug) =  nnz(~real_treatfails_gainedres_toprch(was_reccomended(:,drug)));
   x(1,2,drug) =  nnz(real_treatfails_gainedres_toprch(was_not_reccomended(:,drug)));
   x(2,2,drug) =  nnz(~real_treatfails_gainedres_toprch(was_not_reccomended(:,drug)));
    
    rate_gain_res_rec2(drug) = nnz(real_treatfails_gainedres_toprch(was_reccomended(:,drug)))/length(real_treatfails_gainedres_toprch(was_reccomended(:,drug)));
    total_number_rec2 = length(real_treatfails_gainedres_toprch(was_reccomended(:,drug)));
    
    rate_gain_res_rec2_error(drug) = sqrt((rate_gain_res_rec2(drug)*(1-rate_gain_res_rec2(drug)))/total_number_rec2);
   
    rate_gain_res_notrec(drug) = nnz(real_treatfails_gainedres_toprch(was_not_reccomended(:,drug)))/length(real_treatfails_gainedres_toprch(was_not_reccomended(:,drug)));
    total_number_notrec = length(real_treatfails_gainedres_toprch(was_not_reccomended(:,drug)));
    rate_gain_res_notrec_error(drug) = sqrt((rate_gain_res_notrec(drug)*(1-rate_gain_res_notrec(drug)))/total_number_notrec);
   
    [h1,p1] = fishertest(x(:,:,drug));
    p_fisher(drug) = p1;
       
    n1 = x(1,1,drug); N1 = sum(x(1:2,1,drug));
       n2 = x(1,2,drug); N2 = sum(x(1:2,2,drug));
       % Pooled estimate of proportion
       p0 = (n1+n2) / (N1+N2);
       % Expected counts under H0 (null hypothesis)
       n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];
       chi2stat = sum((observed-expected).^2 ./ expected);
       p(drug) = 1 - chi2cdf(chi2stat,1);

    
    
end


%
figure;
set(gcf,'color','w', 'name','Fig. 4D', 'units','centimeters','Position',[1 1 15 10]);
width = 0.5;
gap = 1.8;
plot([0 7*gap+3], [0 0],'k')

for ii = 1:number_of_drugs
hold on
xpos = ii*gap-0.32;
hd3(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 rate_gain_res_notrec(1,ii)*100 rate_gain_res_notrec(1,ii)*100],'w','EdgeColor','k', 'FaceColor', [0.78 0.13 0.13]); % using MATLAB "patch" to establish the border
errorbar(xpos, rate_gain_res_notrec(1,ii)*100,rate_gain_res_notrec_error(1,ii)*100,'k', 'LineStyle','none')
end

for ii = 1:number_of_drugs
hold on
xpos = ii*gap+0.32;
hd1(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 rate_gain_res_rec2(1,ii)*100 rate_gain_res_rec2(1,ii)*100],'w','EdgeColor','k', 'FaceColor', [0 0.6 0.07]); % using MATLAB "patch" to establish the border
errorbar(xpos, rate_gain_res_rec2(1,ii)*100,rate_gain_res_rec2_error(1,ii)*100,'k', 'LineStyle','none')

end
xlim([0.5 number_of_drugs*gap+1])
ylabel('% of patients who acquired post-treatment resistance ')
xlabel('treatment drug')
set(gca,'XTick',gap:gap:gap*7);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(1:7));
xtickangle(45);

hold on
y = rate_gain_res_notrec*100;
for ii = 1:length(y)
    xpos = ii*gap;
    if p(ii)< 0.05 && p(ii)> 0.005
    plot(xpos, y(ii)*1.1, '*k', 'MarkerSize',4)
    elseif p(ii) <= 0.005 && p(ii)> 0.0005
    plot([xpos-0.1 xpos+0.1] , [y(ii)*1.1 y(ii)*1.1], '*k', 'MarkerSize',4)
    elseif p(ii) <= 0.0005
        plot([xpos-0.2  xpos xpos+0.2] , [ y(ii)*1.1  y(ii)*1.1 y(ii)*1.1], '*k', 'MarkerSize',4)
    end
end

%% number of each antibiotic prescribed in test period for 
% different prescription methods 

number_prescribed_reccomended = zeros(3,size(real_prch,2));
number_prescribed_reccomended(1,:) = sum(real_prch);
number_prescribed_reccomended(2,:) = sum(best_drug_costadjust);
number_prescribed_reccomended(3,:) = sum(best_drug);


figure;
set(gcf,'color','w', 'name','Fig. S14', 'units','centimeters','Position',[1 1 15 10]);
bar(number_prescribed_reccomended')
legend({'physician prescribed drugs' ; 'ML, recommended cost adjusted' ; 'ML, recommended'})
xticklabels(UTI_cases.SMP_Res_drug_names)
xtickangle(45)
ylabel('number of prescriptions/ reccomendations')


%% Look at recommendatoins in antibiotic sensitive cases where there was past reisitance to that antibitoc 
real_prch = UTI_cases.PCR_sameday(withPdata_index, 1:number_of_drugs);
figure;
set(gcf,'color','w', 'name','Fig. S15', 'units','centimeters','Position',[1 1 12 25]);

for drug = 1:7
any_previous_cipr = num_of_previous_R(withPdata_index ,drug)>0 ;
cip_sensitve = ismember(UTI_cases.SMP_Res(withPdata_index,drug),params.sensitive_group);
only_cip = any_previous_cipr & cip_sensitve;


number_prescribed_reccomended = zeros(3,size(real_prch(only_cip,:),2));
number_prescribed_reccomended(1,:) = sum(real_prch(only_cip,:));
number_prescribed_reccomended(2,:) = sum(best_drug_costadjust(only_cip,:));
number_prescribed_reccomended(3,:) = sum(best_drug(only_cip,:));
diff_recc(1,:) = -number_prescribed_reccomended(1,:)+number_prescribed_reccomended(2,:);
diff_recc(2,:) = -number_prescribed_reccomended(1,:)+number_prescribed_reccomended(3,:);
subplot(4,2,drug)
bar(diff_recc');
ylim([-450 450]);
%set(gca,'YScale','log')
xticklabels(UTI_cases.SMP_Res_drug_names);
xtickangle(45);
ylabel({'difference between number of antibiotic'; 'prescribed by physicans and'; 'reccomended by ML models'})
title(UTI_cases.SMP_Res_drug_names{drug})
end

%%
%% was the treatment reccomended or not for differnt tresholds
%bins =0:0.001:0.3;
figure;
set(gcf,'color','w', 'name','Fig. S13', 'units','centimeters','Position',[1 1 40 10]);
thresholds = [90 70 50]; % binarize using 3 different thresholds

for mm = 1:length(thresholds)
clear x
below_thresh = zeros(size(P_nonans));
for drug = 1:number_of_drugs
threshold(drug) = prctile(P_nonans(P_nonans(:,drug)<1,drug),thresholds(mm));
below_thresh(:,drug) = P_nonans(:,drug) < threshold(drug);
below_thresh(P_nonans(:,drug) > 1,drug) = NaN;
end

above_thresh = below_thresh == 0;
below_thresh = below_thresh == 1;

was_reccomended = logical(real_prch) & ismember(real_res, [ 1 2]) & below_thresh;
was_not_reccomended = logical(real_prch) & ismember(real_res, [ 1 2]) & above_thresh;

p = zeros(1,number_of_drugs);
for drug = 1:number_of_drugs
   x(1,1,drug) =  nnz(real_treatfails_gainedres_toprch(was_reccomended(:,drug)));
   x(2,1,drug) =  nnz(~real_treatfails_gainedres_toprch(was_reccomended(:,drug)));
   x(1,2,drug) =  nnz(real_treatfails_gainedres_toprch(was_not_reccomended(:,drug)));
   x(2,2,drug) =  nnz(~real_treatfails_gainedres_toprch(was_not_reccomended(:,drug)));
    
    rate_gain_res_rec2(drug) = nnz(real_treatfails_gainedres_toprch(was_reccomended(:,drug)))/length(real_treatfails_gainedres_toprch(was_reccomended(:,drug)));
    total_number_rec2 = length(real_treatfails_gainedres_toprch(was_reccomended(:,drug)));
    
    rate_gain_res_rec2_error(drug) = sqrt((rate_gain_res_rec2(drug)*(1-rate_gain_res_rec2(drug)))/total_number_rec2);
   
    rate_gain_res_notrec(drug) = nnz(real_treatfails_gainedres_toprch(was_not_reccomended(:,drug)))/length(real_treatfails_gainedres_toprch(was_not_reccomended(:,drug)));
    total_number_notrec = length(real_treatfails_gainedres_toprch(was_not_reccomended(:,drug)));
    rate_gain_res_notrec_error(drug) = sqrt((rate_gain_res_notrec(drug)*(1-rate_gain_res_notrec(drug)))/total_number_notrec);
   
    [h1,p1] = fishertest(x(:,:,drug));
    p_fisher(drug) = p1;
       
    n1 = x(1,1,drug); N1 = sum(x(1:2,1,drug));
       n2 = x(1,2,drug); N2 = sum(x(1:2,2,drug));
       % Pooled estimate of proportion
       p0 = (n1+n2) / (N1+N2);
       % Expected counts under H0 (null hypothesis)
       n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];
       chi2stat = sum((observed-expected).^2 ./ expected);
       p(drug) = 1 - chi2cdf(chi2stat,1);

    
    
end


%
subplot(1,3,mm)
width = 0.5;
gap = 1.8;
plot([0 7*gap+3], [0 0],'k')

for ii = 1:number_of_drugs
hold on
xpos = ii*gap-0.32;
hd3(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 rate_gain_res_notrec(1,ii)*100 rate_gain_res_notrec(1,ii)*100],'w','EdgeColor','k', 'FaceColor', [0.78 0.13 0.13]); % using MATLAB "patch" to establish the border
errorbar(xpos, rate_gain_res_notrec(1,ii)*100,rate_gain_res_notrec_error(1,ii)*100,'k', 'LineStyle','none')

%errorbar(xpos,(fail_split_noprch(1,ii)+fail_split_noprch(2,ii)),ratio_noprch_error(ii), 'k', 'LineStyle','none')

end

for ii = 1:number_of_drugs
hold on
xpos = ii*gap+0.32;
hd1(ii) = patch([xpos-width/2 xpos+width/2 xpos+width/2 xpos-width/2],[0  0 rate_gain_res_rec2(1,ii)*100 rate_gain_res_rec2(1,ii)*100],'w','EdgeColor','k', 'FaceColor', [0 0.6 0.07]); % using MATLAB "patch" to establish the border
errorbar(xpos, rate_gain_res_rec2(1,ii)*100,rate_gain_res_rec2_error(1,ii)*100,'k', 'LineStyle','none')

end
xlim([0.5 number_of_drugs*gap+1])
title(['threshold' num2str(100-thresholds(mm))]);
ylabel('% of patients who acquired post-treatment resistance ')
xlabel('treatment drug')
set(gca,'XTick',gap:gap:gap*7);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(1:7));
xtickangle(45);

hold on
y = rate_gain_res_notrec*100;
for ii = 1:length(y)
    xpos = ii*gap;
    if p(ii)< 0.05 && p(ii)> 0.005
    plot(xpos, y(ii)*1.1, '*k', 'MarkerSize',4)
    elseif p(ii) <= 0.005 && p(ii)> 0.0005
    plot([xpos-0.1 xpos+0.1] , [y(ii)*1.1 y(ii)*1.1], '*k', 'MarkerSize',4)
    elseif p(ii) <= 0.0005
        plot([xpos-0.2  xpos xpos+0.2] , [ y(ii)*1.1  y(ii)*1.1 y(ii)*1.1], '*k', 'MarkerSize',4)
    end
end

end

%% Do regression for any failure (S-S & S-R) and look at the predicted risk for the 
% antibiotic recomendations trained to miniize StoR failures
sen_res = [ 1 2 3];
clear coef

for drug = 1:7

SMP_to_use =  find((UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & ~test_data & UTI_cases.hasdiag);
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);

X = [X_demog_to_use(drug_prch_temp==1,:), num_previous_R(drug_prch_temp==1), num_previous_S(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
[glm.B,glm.dev,glm.stats] = glmfit(X,Y,'binomial','logit') ;
coef(:,drug) = glm.B;

end
%% test data:
clear P
P = nan(height(UTI_cases.Demog),number_of_drugs);
treatfails_gainedres_toprch = zeros(height(UTI_cases.Demog),number_of_drugs);
%patient was correctly prescribed a drug from the used drugs
was_correctly_prescribed = sum(UTI_cases.PCR_sameday(:,1:number_of_drugs) & ismember(UTI_cases.SMP_Res(:,1:number_of_drugs), params.sensitive_group),2) > 0;

for drug = 1:number_of_drugs

SMP_to_use =  find(was_correctly_prescribed  & test_data & UTI_cases.hasdiag);
num_previous_R = num_of_previous_R(SMP_to_use ,drug);
num_previous_S = num_of_previous_S(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use, num_previous_R num_previous_S];
P(SMP_to_use, drug) = glmval(coef(:,drug),X,'logit') ;

%%
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
treatfails = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
treatfails_gainedres_toprch(SMP_to_use, drug) = treatfails ;%;& drug_prch_temp;

%% reassign any mismatched (or unknown Res) probabilities to 1
index = find(ismember(UTI_cases.SMP_Res(SMP_to_use,drug), [0 3]));
P(SMP_to_use(index), drug)  = 1.1; % mismatched treatments are not used 

end

%
withPdata_index_SS = find(sum(isnan(P),2)== 0);
P_nonans_SS = P(withPdata_index_SS,:);

%using the best_drug etc from the gain of resistance S-R model  
%on the P values from the S-any model
best_cost_SS_boot = bootstrp(num_bootstraps,@mean,P_nonans_SS(best_drug_costadjust ==1));
best_SS_boot = bootstrp(num_bootstraps,@mean,P_nonans_SS(best_drug ==1));
real_SS_boot = bootstrp(num_bootstraps,@mean,P_nonans_SS(real_prch ==1));


figure;
set(gcf,'color','w', 'name','Fig. S17 Overall risk of recurrence', 'units','centimeters','Position',[1 1 10 10]);
b1= bar([mean(real_SS_boot) mean(best_SS_boot) mean(best_cost_SS_boot) ]);
hold on
errorbar(b1(1).XData,b1(1).YData,[std(real_SS_boot) std(best_SS_boot) std(best_cost_SS_boot)  ],'k', 'LineStyle','none')
xticklabels({'All rec real','All rec best','All rec cost'})
xtickangle( 45 )
ylabel('predicted probability of recurrence')
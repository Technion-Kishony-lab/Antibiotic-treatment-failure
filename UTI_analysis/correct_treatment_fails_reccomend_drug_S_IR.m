function [] = correct_treatment_fails_reccomend_drug_S_IR(UTI_cases,params)
% function takes the UTI_case structure and optional params and treains and tests 
% regressoin model to predict risk of gain of resistance for all susceptibility-matched treated cases cases 
% with the intermediate grouped with resistant
X_age = UTI_cases.Demog.Age;
X_age(:,6) = []; % reference age
X_demog = [X_age, UTI_cases.Demog.Gender, UTI_cases.Demog.Preg UTI_cases.Demog.any_prev_cath];
number_of_drugs = 7; % not including orflox as not measured in test period

%% seperate training data and test data
test_date_start = datenum('2018-05-01') ;
test_date_end = datenum('2019-07-01') ;
test_date_range = test_date_start:test_date_end;
test_data = ismember(UTI_cases.SamplingDate, test_date_range);
fprintf('Testing using %.2f percent of the data\n' ,nnz(test_data)/length(test_data)*100)

%% Regression for StoR failure given any previous resistance to drug
clear coef
num_of_previous_R = UTI_cases.num_of_previous_R;
num_of_previous_S = UTI_cases.num_of_previous_S;
sen_res = [ 2 3]; % intermediate grouped with resistant

for drug = 1:number_of_drugs
    
SMP_to_use =  find( ismember(UTI_cases.SMP_Res(:,drug), 1)  & ~test_data & UTI_cases.hasdiag);
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
was_correctly_prescribed = sum(UTI_cases.PCR_sameday(:,1:number_of_drugs) & ismember(UTI_cases.SMP_Res(:,1:number_of_drugs), [1]),2) > 0;

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
treatfails_gainedres_toprch(SMP_to_use, drug) = treatfails & drug_prch_temp;

%% reassign any mismatched (or unknown Res) probabilities to 1
index = find(ismember(UTI_cases.SMP_Res(SMP_to_use,drug), [0 3]));
P(SMP_to_use(index), drug)  = 1.1;

end

treatfails_gainedres_toprch = logical(sum(treatfails_gainedres_toprch,2));


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

%% was the treatment reccomended or not
clear x
below_thresh = zeros(size(P_nonans));
for drug = 1:number_of_drugs
threshold(drug) = prctile(P_nonans(P_nonans(:,drug)<1,drug),85);
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

subplot(1,2,2)
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
%ylim([0 23])
ylabel('% of patients who acquired post-treatment resistance ')
xlabel('treatment drug')
% set(gca,'XTick',2:2:16);
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

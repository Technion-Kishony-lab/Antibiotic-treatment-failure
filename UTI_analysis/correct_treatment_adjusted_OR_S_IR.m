function [] = correct_treatment_adjusted_OR_S_IR(UTI_cases,params)
% function takes the UTI_case structure and optional params and performs
% regresion for each antibitoic to calulate demographics adjusted odds
% ratio of risk of early recurrence, for cases susceptibitlyi-matched treated cases which gained resistance
% and remained sensitive. Intermediate now grouped with resistant

X_age = UTI_cases.Demog.Age;
X_age(:,6) = []; % refence
X_demog = [X_age, UTI_cases.Demog.Gender, UTI_cases.Demog.Preg UTI_cases.Demog.any_prev_cath ];
number_of_drugs = 7; %not including ofloxacin for consistency with drug recommendations 
min_1 = 10;

%% OR specifically StoR failure given any previous resistance to drug
clear y_SSR p_SSR ciu_SSR cil_SSR se_SSR 

sen_res = [2 3]; % intermediate grouped with resistant

seall = zeros(size(X_demog,2) +2, number_of_drugs);
yall = zeros(size(X_demog,2) +2, number_of_drugs);
pall = zeros(size(X_demog,2) +2, number_of_drugs);

for drug = 1:number_of_drugs

SMP_to_use =  find(UTI_cases.any_SRmeasurement(:,drug)>0 & ismember(UTI_cases.SMP_Res(:,drug), 1) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag);
any_previous_R = UTI_cases.num_of_previous_R(SMP_to_use ,drug)>0 ;
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), any_previous_R(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
c = my_fit(X, Y, min_1);
seall(:,drug) = c.se;
yall(:,drug) = c.coef;
pall(:, drug) = c.p;
end

y_SSR(2:2:number_of_drugs*2) = exp(yall(end,:));
p_SSR(2:2:number_of_drugs*2) = pall(end,:);
se_SSR(2:2:number_of_drugs*2) = seall(end,:);
cil_SSR(2:2:number_of_drugs*2) = exp(yall(end,:) -  seall(end,:));
ciu_SSR(2:2:number_of_drugs*2) = exp(yall(end,:) +  seall(end,:));

x1 = yall(end,:);
se1 = seall(end,:);


%% OR specifically StoR failure given any previous resistance to drug
sen_res = [1]; % intermediate not grouped with sensitive

seall = zeros(size(X_demog,2) +2, number_of_drugs);
yall = zeros(size(X_demog,2) +2, number_of_drugs);
pall = zeros(size(X_demog,2) +2, number_of_drugs);

for drug = 1:number_of_drugs
SMP_to_use =  find(UTI_cases.any_SRmeasurement(:,drug)>0 & ismember(UTI_cases.SMP_Res(:,drug),1) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag);
any_previous_R = UTI_cases.num_of_previous_R(SMP_to_use ,drug)>0 ;
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), any_previous_R(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
c = my_fit(X, Y,min_1);
seall(:,drug) = c.se;
yall(:,drug) = c.coef;
pall(:, drug) = c.p;

end

y_SSR(1:2:number_of_drugs*2-1) = exp(yall(end,:));
p_SSR(1:2:number_of_drugs*2-1) = pall(end,:);
se_SSR(1:2:number_of_drugs*2-1) = seall(end,:);

cil_SSR(1:2:number_of_drugs*2-1) = exp(yall(end,:) -  seall(end,:));
ciu_SSR(1:2:number_of_drugs*2-1) = exp(yall(end,:) +  seall(end,:));

b1 = bar(y_SSR(1:number_of_drugs*2), 'FaceColor','flat',  'BarWidth', 0.8);
hold on

b1(1).CData(1:2:number_of_drugs*2-1,:) = ones(number_of_drugs,1)*params.SS_color;
b1.CData(2:2:number_of_drugs*2,:) = ones(number_of_drugs,1)*params.SR_color;
errneg = y_SSR - cil_SSR;
errpos =  ciu_SSR -y_SSR;
errorbar(b1.XData,y_SSR,errneg,errpos,'k', 'LineStyle','none');

hold on
for ii = 1:length(y_SSR)
    if p_SSR(ii)< 0.01 && p_SSR(ii)> 0.001
    plot(b1.XData(ii), y_SSR(ii)+1, '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.001 && p_SSR(ii)> 0.00001
    plot([b1.XData(ii)-0.1 b1.XData(ii)+0.1] , [y_SSR(ii)+1 y_SSR(ii)+1], '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.00001
        plot([b1.XData(ii)-0.2  b1.XData(ii) b1.XData(ii)+0.2] , [ y_SSR(ii)+1  y_SSR(ii)+1 y_SSR(ii)+1], '*k', 'MarkerSize',4)
    end
end

xticks((1:2:number_of_drugs*2)+0.5);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(1:number_of_drugs));
xtickangle(45);

ylabel({'Adjusted odds ratio of early UTI recurrence';...
    'given any past resistant sample'});

text(0,6.5,'{\it * P<0.01}');
text(0,6,'{\it** P<0.001}');
text(0,5.5,'{\it*** P<0.00001}');

set(gca, 'YScale', 'log');
ylim([1 7])

end
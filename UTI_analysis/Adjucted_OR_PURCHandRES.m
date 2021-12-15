function [] = Adjucted_OR_PURCHandRES(UTI_cases,params)
% function takes the UTI_case structure and optional params and performs
% regresion for each antibitoic to calulate demographics adjusted odds
% ratio of risk of early recurrence, for cases susceptibitlyi-matched
% treated cases including both past infection susceptibiltiy and past
% antibiotic purchases

%Adjusted for following demographics:
X_age = UTI_cases.Demog.Age;
X_age(:,6) = []; %reference age
X_demog = [X_age, UTI_cases.Demog.Gender, UTI_cases.Demog.Preg UTI_cases.Demog.any_prev_cath ];
number_of_drugs = 7; %not including ofloxacin for consistency with drug recommendatoins 
num_prev_purch = UTI_cases.num_prev_purch; % also now include num prev purch of treated drug

min_1 = 10;
%% OR specifically StoR failure given any previous resistance to drug
clear y_SSR p_SSR ciu_SSR cil_SSR se_SSR 

sen_res = params.resistant_group; % gained res
seall = zeros(size(X_demog,2) +3, number_of_drugs);
yall = zeros(size(X_demog,2) +3, number_of_drugs);
pall = zeros(size(X_demog,2) +3, number_of_drugs);

for drug = 1:number_of_drugs
SMP_to_use =  find(UTI_cases.any_SRmeasurement(:,drug)>0 & (UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag);
any_previous_purch = num_prev_purch(SMP_to_use ,drug)>0 ;
any_previous_R = UTI_cases.num_of_previous_R(SMP_to_use ,drug)>0 ;
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), any_previous_R(drug_prch_temp==1) any_previous_purch(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
c = my_fit(X, Y, min_1);
seall(:,drug) = c.se;
yall(:,drug) = c.coef;
pall(:, drug) = c.p;
end

y_SSR(2:2:number_of_drugs*2) = exp(yall(end-1,:));
p_SSR(2:2:number_of_drugs*2) = pall(end-1,:);
se_SSR(2:2:number_of_drugs*2) = seall(end-1,:);
cil_SSR(2:2:number_of_drugs*2) = exp(yall(end-1,:) -  seall(end-1,:));
ciu_SSR(2:2:number_of_drugs*2) = exp(yall(end-1,:) +  seall(end-1,:));

x1 = yall(end,:)
se1 = seall(end,:)

y_SSR_purch(2:2:number_of_drugs*2) = exp(yall(end,:));
p_SSR_purch(2:2:number_of_drugs*2) = pall(end,:);
se_SSR_purch(2:2:number_of_drugs*2) = seall(end,:);
cil_SSR_purch(2:2:number_of_drugs*2) = exp(yall(end,:) -  seall(end,:));
ciu_SSR_purch(2:2:number_of_drugs*2) = exp(yall(end,:) +  seall(end,:));

x1_purch = yall(end,:);
se1_purch = seall(end,:);


%% remained Sensitive
sen_res = params.sensitive_group;
seall = zeros(size(X_demog,2) +3, number_of_drugs);
yall = zeros(size(X_demog,2) +3, number_of_drugs);
pall = zeros(size(X_demog,2) +3, number_of_drugs);

for drug = 1:number_of_drugs
SMP_to_use =  find(UTI_cases.any_SRmeasurement(:,drug)>0 & (UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag);
any_previous_purch = num_prev_purch(SMP_to_use ,drug)>0 ;
any_previous_R = UTI_cases.num_of_previous_R(SMP_to_use ,drug)>0 ;
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), any_previous_R(drug_prch_temp==1) any_previous_purch(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
c = my_fit(X, Y,min_1);
seall(:,drug) = c.se;
yall(:,drug) = c.coef;
pall(:, drug) = c.p;
end
%
y_SSR(1:2:number_of_drugs*2-1) = exp(yall(end-1,:));
p_SSR(1:2:number_of_drugs*2-1) = pall(end-1,:);
se_SSR(1:2:number_of_drugs*2-1) = seall(end-1,:);
cil_SSR(1:2:number_of_drugs*2-1) = exp(yall(end-1,:) -  seall(end-1,:));
ciu_SSR(1:2:number_of_drugs*2-1) = exp(yall(end-1,:) +  seall(end-1,:));

y_SSR_purch(1:2:number_of_drugs*2) = exp(yall(end,:));
p_SSR_purch(1:2:number_of_drugs*2) = pall(end,:);
se_SSR_purch(1:2:number_of_drugs*2) = seall(end,:);
cil_SSR_purch(1:2:number_of_drugs*2) = exp(yall(end,:) -  seall(end,:));
ciu_SSR_purch(1:2:number_of_drugs*2) = exp(yall(end,:) +  seall(end,:));

%%
subplot(2,2,1)
b1 = bar(y_SSR(1:number_of_drugs*2), 'FaceColor','flat',  'BarWidth', 0.8);
hold on

b1(1).CData(1:2:number_of_drugs*2-1,:) = ones(number_of_drugs,1)*params.SS_color;
b1.CData(2:2:number_of_drugs*2,:) = ones(number_of_drugs,1)*params.SR_color;
%plot(b1.XData,cil_SSR,'kx');
%plot(b1.XData,ciu_SSR,'kx');
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
% c.se(c.p> 0.05) = 0;
% hold on
% errorbar(c.coef(2:end),c.se(2:end),'.')
xticks((1:2:number_of_drugs*2)+0.5);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(1:number_of_drugs));
xtickangle(45);

ylabel({'Adjusted odds ratio of early recurrence'; ...
    'given past resitant sample'}, 'FontSize' , 8);

text(0,6.5,'{\it * P<0.01}');
text(0,6,'{\it** P<0.001}');
text(0,5.5,'{\it*** P<0.00001}');

set(gca, 'YScale', 'log');
ylim([1 7])
title('past resistant isolates')


%%
subplot(2,2,2)
b1 = bar(y_SSR_purch(1:number_of_drugs*2), 'FaceColor','flat',  'BarWidth', 0.8);
hold on

b1(1).CData(1:2:number_of_drugs*2-1,:) = ones(number_of_drugs,1)*params.SS_color;
b1.CData(2:2:number_of_drugs*2,:) = ones(number_of_drugs,1)*params.SR_color;
%plot(b1.XData,cil_SSR,'kx');
%plot(b1.XData,ciu_SSR,'kx');
errneg = y_SSR_purch - cil_SSR_purch;
errpos =  ciu_SSR_purch -y_SSR_purch;

errorbar(b1.XData,y_SSR_purch,errneg,errpos,'k', 'LineStyle','none');

hold on
for ii = 1:length(y_SSR_purch)
    if p_SSR_purch(ii)< 0.01 && p_SSR_purch(ii)> 0.001
    plot(b1.XData(ii), y_SSR_purch(ii)+1, '*k', 'MarkerSize',4)
    elseif p_SSR_purch(ii) <= 0.001 && p_SSR(ii)> 0.00001
    plot([b1.XData(ii)-0.1 b1.XData(ii)+0.1] , [y_SSR_purch(ii)+1 y_SSR_purch(ii)+1], '*k', 'MarkerSize',4)
    elseif p_SSR_purch(ii) <= 0.00001
        plot([b1.XData(ii)-0.2  b1.XData(ii) b1.XData(ii)+0.2] , [ y_SSR_purch(ii)+1  y_SSR_purch(ii)+1 y_SSR_purch(ii)+1], '*k', 'MarkerSize',4)
    end
end
% c.se(c.p> 0.05) = 0;
% hold on
% errorbar(c.coef(2:end),c.se(2:end),'.')
xticks((1:2:number_of_drugs*2)+0.5);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(1:number_of_drugs));
xtickangle(45);

ylabel({'Adjusted odds ratio of early recurrence'; ...
    'given a past antibitoic purchase'}, 'FontSize' , 8);

text(0,6.5,'{\it * P<0.01}');
text(0,6,'{\it** P<0.001}');
text(0,5.5,'{\it*** P<0.00001}');

set(gca, 'YScale', 'log');

ylim([1 7])
title('past purchases')
end
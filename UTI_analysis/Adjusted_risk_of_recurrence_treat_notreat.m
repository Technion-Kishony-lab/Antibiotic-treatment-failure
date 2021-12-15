function [] = Adjusted_risk_of_recurrence_treat_notreat(UTI_cases,params)
% function takes the UTI_case structure and optional params and performs
% regresion to calulate demographics adjusted odds ratio of risk of early recurrence
% given treatment with specific antiboitic compared to untreated

number_of_drugs = params.number_drugs;
X_age = UTI_cases.Demog.Age;
X_age(:,6) = []; %reference age
X_demog = [X_age, UTI_cases.Demog.Gender, UTI_cases.Demog.Preg ];


%% OR of StoR failure given treated vs untreated
clear y_SSR p_SSR ciu_SSR cil_SSR se_SSR 
%exclude treatment fails less than 5 days

sen_res = params.resistant_group;

seall_gained = zeros(size(X_demog,2) +2, number_of_drugs);
yall_gained = zeros(size(X_demog,2) +2, number_of_drugs);
pall_gained = zeros(size(X_demog,2) +2, number_of_drugs);

for drug = 1:number_of_drugs

SMP_to_use =  find((UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag & (UTI_cases.PCR_sameday(:,drug)==1 | UTI_cases.PCR_sameday(:,10)==1));
treated = UTI_cases.PCR_sameday(SMP_to_use ,drug);
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug) | UTI_cases.PCR_sameday(SMP_to_use ,10);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), treated(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
c = my_fit(X, Y, 50);
seall_gained(:,drug) = c.se;
yall_gained(:,drug) = c.coef;
pall_gained(:, drug) = c.p;

end
%

yall_gained(isnan(yall_gained)) = 0;
y_SSR(1:number_of_drugs) = exp(yall_gained(end,:));
p_SSR(1:number_of_drugs) = pall_gained(end,:);
se_SSR(1:number_of_drugs) = seall_gained(end,:);
cil_SSR(1:number_of_drugs) = exp(yall_gained(end,:) -  seall_gained(end,:));
ciu_SSR(1:number_of_drugs) = exp(yall_gained(end,:) +  seall_gained(end,:));

subplot(3,1,1)
b1 = bar(y_SSR(1:number_of_drugs), 'FaceColor','flat',  'BarWidth', 0.8);
hold on

b1(1).CData(1:number_of_drugs,:) = ones(number_of_drugs,1)*params.SR_color;
errneg = y_SSR - cil_SSR;
errpos =  ciu_SSR -y_SSR;
errorbar(b1.XData,y_SSR,errneg,errpos,'k', 'LineStyle','none');

hold on
above_gap = 0.5;
for ii = 1:length(y_SSR)
    if p_SSR(ii)< 0.01 && p_SSR(ii) > 0.001
    plot(b1.XData(ii), y_SSR(ii)+above_gap, '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.001 && p_SSR(ii)> 0.0001
    plot([b1.XData(ii)-0.1 b1.XData(ii)+0.1] , [y_SSR(ii)+above_gap y_SSR(ii)+above_gap], '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.0001
        plot([b1.XData(ii)-0.2  b1.XData(ii) b1.XData(ii)+0.2] , [ y_SSR(ii)+above_gap  y_SSR(ii)+above_gap y_SSR(ii)+above_gap], '*k', 'MarkerSize',4)
    end
end

xticks(1:number_of_drugs);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(1:number_of_drugs));
xtickangle(45);
ylabel({'Adjusted odds ratio of infections sensitive';...
    'to the focal antibiotic gaining resistance when';...
    ' either treated with the focal antibiotic or untreated'}, 'FontSize', 8);
xlabel('focal antibiotic')
set(gca, 'YScale', 'log');
ylim([1 6])

%% OR of any failure given treated vs untreated for intially sensitive cases
seall = zeros(size(X_demog,2) +2, number_of_drugs);
yall = zeros(size(X_demog,2) +2, number_of_drugs);
pall = zeros(size(X_demog,2) +2, number_of_drugs);

subplot(3,1,2)
sen_res = [1 2 3]; % any recurrence


for drug = 1:number_of_drugs

SMP_to_use =  find((UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag & (UTI_cases.PCR_sameday(:,drug)==1 | UTI_cases.PCR_sameday(:,10)==1));
treated = UTI_cases.PCR_sameday(SMP_to_use ,drug);
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug) | UTI_cases.PCR_sameday(SMP_to_use ,10);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), treated(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
c = my_fit(X, Y, 50);
seall(:,drug) = c.se;
yall(:,drug) = c.coef;
pall(:, drug) = c.p;
end

yall(isnan(yall)) = 0;
y_SSR(1:number_of_drugs) = exp(yall(end,:));
p_SSR(1:number_of_drugs) = pall(end,:);
se_SSR(1:number_of_drugs) = seall(end,:);
cil_SSR(1:number_of_drugs) = exp(yall(end,:) -  seall(end,:));
ciu_SSR(1:number_of_drugs) = exp(yall(end,:) +  seall(end,:));
b1 = bar(y_SSR(1:number_of_drugs), 'FaceColor','flat',  'BarWidth', 0.8);
hold on
b1(1).CData(1:number_of_drugs,:) = ones(number_of_drugs,1)*[0.4 0.4 0.4];
errneg = y_SSR - cil_SSR;
errpos =  ciu_SSR -y_SSR;
errorbar(b1.XData,y_SSR,errneg,errpos,'k', 'LineStyle','none');

hold on
above_gap = 0.05;
for ii = 1:length(y_SSR)
    if p_SSR(ii)< 0.01 && p_SSR(ii)> 0.01
    plot(b1.XData(ii), y_SSR(ii)+above_gap, '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.01 && p_SSR(ii)> 0.0001
    plot([b1.XData(ii)-0.1 b1.XData(ii)+0.1] , [y_SSR(ii)+above_gap y_SSR(ii)+above_gap], '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.0001
        plot([b1.XData(ii)-0.2  b1.XData(ii) b1.XData(ii)+0.2] , [ y_SSR(ii)+above_gap  y_SSR(ii)+above_gap y_SSR(ii)+above_gap], '*k', 'MarkerSize',4)
    end
end

xticks(1:number_of_drugs);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(1:number_of_drugs));
xtickangle(45);
ylabel({'Adjusted odds ratio of infections sensitive';...
    'to the focal antibiotic resulting in early recurrence when';...
    ' either treated with the focal antibiotic or untreated'}, 'FontSize', 8);
xlabel('focal antibiotic')
text(9,0.15,'{\it * P<0.01}');
text(9,0.14,'{\it** P<0.0001}');
text(9,0.13,'{\it*** P<0.00001}');
set(gca, 'YScale', 'log');
ylim([0.2 1.4])

%% OR of any failure given treated vs untreated for intially resistant cases
seall = zeros(size(X_demog,2) +2, number_of_drugs);
yall = zeros(size(X_demog,2) +2, number_of_drugs);
pall = zeros(size(X_demog,2) +2, number_of_drugs);
subplot(3,1,3)
sen_res = [1 2 3]; % any recurrence

for drug = 1:number_of_drugs

SMP_to_use =  find(UTI_cases.SMP_Res(:,drug)== 3 & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag & (UTI_cases.PCR_sameday(:,drug)==1 | UTI_cases.PCR_sameday(:,10)==1));
treated = UTI_cases.PCR_sameday(SMP_to_use ,drug);
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug) | UTI_cases.PCR_sameday(SMP_to_use ,10);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), treated(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
c = my_fit(X, Y, 50);
seall(:,drug) = c.se;
yall(:,drug) = c.coef;
pall(:, drug) = c.p;
end
yall(isnan(yall)) = 0;
y_SSR(1:number_of_drugs) = exp(yall(end,:));
p_SSR(1:number_of_drugs) = pall(end,:);
se_SSR(1:number_of_drugs) = seall(end,:);
cil_SSR(1:number_of_drugs) = exp(yall(end,:) -  seall(end,:));
ciu_SSR(1:number_of_drugs) = exp(yall(end,:) +  seall(end,:));
b1 = bar(y_SSR(1:number_of_drugs), 'FaceColor','flat',  'BarWidth', 0.8);
hold on
b1(1).CData(1:number_of_drugs,:) = ones(number_of_drugs,1)*[0.75 0.75 0.75];
errneg = y_SSR - cil_SSR;
errpos =  ciu_SSR -y_SSR;
errorbar(b1.XData,y_SSR,errneg,errpos,'k', 'LineStyle','none');
hold on
above_gap = 0.05;
for ii = 1:length(y_SSR)
    if p_SSR(ii)< 0.01 && p_SSR(ii)> 0.01
    plot(b1.XData(ii), y_SSR(ii)+above_gap, '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.01 && p_SSR(ii)> 0.0001
    plot([b1.XData(ii)-0.1 b1.XData(ii)+0.1] , [y_SSR(ii)+above_gap y_SSR(ii)+above_gap], '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.0001
        plot([b1.XData(ii)-0.2  b1.XData(ii) b1.XData(ii)+0.2] , [ y_SSR(ii)+above_gap  y_SSR(ii)+above_gap y_SSR(ii)+above_gap], '*k', 'MarkerSize',4)
    end
end
xticks(1:number_of_drugs);
set(gca,'xticklabel',UTI_cases.SMP_Res_drug_names(1:number_of_drugs));
xtickangle(45);
ylabel({'Adjusted odds ratio of infections resistant';...
    'to the focal antibiotic resulting in early recurrence when';...
    ' either treated with the focal antibiotic or untreated'}, 'FontSize', 8);
xlabel('focal antibiotic')
text(9,0.15,'{\it * P<0.01}');
text(9,0.14,'{\it** P<0.0001}');
text(9,0.13,'{\it*** P<0.00001}');
set(gca, 'YScale', 'log');
ylim([0.2 1.4])

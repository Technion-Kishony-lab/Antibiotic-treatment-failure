function [] = correct_treatment_adjusted_OR(UTI_cases,params)
% function takes the UTI_case structure and optional params and performs
% regresion for each antibitoic to calulate demographics adjusted odds
% ratio of risk of early recurrence, for cases susceptibitlyi-matched treated cases which gained resistance
% and remained sensitive. Fuction writes table of regressoin coefficients
% to the Tables directory.

%Adjusted for following demographics:
X_age = UTI_cases.Demog.Age;
X_age(:,6) = []; %reference age
X_demog = [X_age, UTI_cases.Demog.Gender, UTI_cases.Demog.Preg UTI_cases.Demog.any_prev_cath ];
features = {'Age 0-9', 'Age 10-19', 'Age 20-29', 'Age 30-39', 'Age 40-49', 'Age 60-69', 'Age 70-79', 'Age 80-89', 'Age 90-100', 'Gender', 'Pregnancy', 'Catheter', 'Prev Resistance'};

number_of_drugs = 7; %not including ofloxacin for consistency with drug recommendatoins 
%% OR specifically StoR failure given any previous resistance to drug
clear y_SSR p_SSR ciu_SSR cil_SSR se_SSR 

sen_res = params.resistant_group;

seall = zeros(size(X_demog,2) +2, number_of_drugs);
yall = zeros(size(X_demog,2) +2, number_of_drugs);
pall = zeros(size(X_demog,2) +2, number_of_drugs);
min_1 = 10;

for drug = 1:number_of_drugs
%use cases sentivie to antibiotic      
SMP_to_use =  find(UTI_cases.any_SRmeasurement(:,drug)>0 & (UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag);
any_previous_R = UTI_cases.num_of_previous_R(SMP_to_use ,drug)>0 ;
% gain of resistance early recurrences
treatfailure = UTI_cases.treatfailure(SMP_to_use) ==1 & ismember(UTI_cases.next_res(SMP_to_use, drug),sen_res) ;
%treated with antibiotic
drug_prch_temp = UTI_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), any_previous_R(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
%do regression
c = my_fit(X, Y, min_1);
seall(:,drug) = c.se;
yall(:,drug) = c.coef;
pall(:, drug) = c.p;
end

y_SSR(2:2:number_of_drugs*2) = exp(yall(end,:)); %OR
p_SSR(2:2:number_of_drugs*2) = pall(end,:);
se_SSR(2:2:number_of_drugs*2) = seall(end,:);
cil_SSR(2:2:number_of_drugs*2) = exp(yall(end,:) -  seall(end,:));
ciu_SSR(2:2:number_of_drugs*2) = exp(yall(end,:) +  seall(end,:));

x1 = yall(end,:);
se1 = seall(end,:);

%% Write TableS4
features = pad(features);
antibiotic_names = pad(UTI_cases.SMP_Res_drug_names);
yall = yall(2:end,:);
pall = pall(2:end,:);
seall= seall(2:end,:);
CIL = yall - seall;
CIU = yall + seall;

pstrs = cell(size(pall)) ;
for i = 1:size(pall,1)
    for j = 1:size(pall,2)
        pstrs{i,j} = '' ;
        if pall(i,j)<0.01  , pstrs{i,j} = [pstrs{i,j}, '*']; end
        if pall(i,j)<0.001 , pstrs{i,j} = [pstrs{i,j}, '*']; end
        if pall(i,j)<0.0001, pstrs{i,j} = [pstrs{i,j}, '*']; end
    end
end

drug_order = [1 2; 3 4; 5 6; 7 7];
fid = fopen('Tables/Table_S7_AdgustedOR_SR.txt', 'w');

for jj = 1:size(drug_order,1)
    
fprintf(fid,'%-12s','');
for drug = drug_order(jj,:)
    fprintf(fid,antibiotic_names{drug});
    fprintf(fid,'%-12s','');
end
fprintf(fid,'\n') ;
for ii = 1:size(yall,1)
    fprintf(fid,[features{ii} '  ']);
for drug = drug_order(jj,:)
    coefs_tmp = [yall(ii, drug) CIL(ii, drug) CIU(ii, drug) pstrs(ii, drug)];
    fprintf(fid,'%5.2f [%5.2f,%5.2f]%-3s  ',coefs_tmp{:}) ;
    fprintf(fid,'%-7s','');
end
fprintf(fid,'\n') ;
end
fprintf(fid,'\n') ;
end
fprintf(fid,'* P<0.01    ** P<0.0001    *** P<0.000001\n\n\n');
fclose(fid) ;

%% OR specifically StoR failure given any previous resistance to drug
sen_res = params.sensitive_group;

seall = zeros(size(X_demog,2) +2, number_of_drugs);
yall = zeros(size(X_demog,2) +2, number_of_drugs);
pall = zeros(size(X_demog,2) +2, number_of_drugs);

for drug = 1:number_of_drugs
SMP_to_use =  find(UTI_cases.any_SRmeasurement(:,drug)>0 & (UTI_cases.SMP_Res(:,drug)== 1 | UTI_cases.SMP_Res(:,drug)== 2) & UTI_cases.next_res(:, drug) ~=0  & UTI_cases.hasdiag);
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


%% TableS8

yall = yall(2:end,:);
pall = pall(2:end,:);
seall= seall(2:end,:);
CIL = yall - seall;
CIU = yall + seall;

pstrs = cell(size(pall)) ;
for i = 1:size(pall,1)
    for j = 1:size(pall,2)
        pstrs{i,j} = '' ;
        if pall(i,j)<0.01  , pstrs{i,j} = [pstrs{i,j}, '*']; end
        if pall(i,j)<0.001 , pstrs{i,j} = [pstrs{i,j}, '*']; end
        if pall(i,j)<0.0001, pstrs{i,j} = [pstrs{i,j}, '*']; end
    end
end

drug_order = [1 2; 3 4; 5 6; 7 7];
fid = fopen('Tables/Table_S8_AdgustedOR_SS.txt', 'w');

for jj = 1:size(drug_order,1)
    
fprintf(fid,'%-12s','');
for drug = drug_order(jj,:)
    fprintf(fid,antibiotic_names{drug});
    fprintf(fid,'%-12s','');
end
fprintf(fid,'\n') ;
for ii = 1:size(yall,1)
    fprintf(fid,[features{ii} '  ']);
for drug = drug_order(jj,:)
    coefs_tmp = [yall(ii, drug) CIL(ii, drug) CIU(ii, drug) pstrs(ii, drug)];
    fprintf(fid,'%5.2f [%5.2f,%5.2f]%-3s  ',coefs_tmp{:}) ;
    fprintf(fid,'%-7s','');
end
fprintf(fid,'\n') ;
end
fprintf(fid,'\n') ;
end
fprintf(fid,'* P<0.01    ** P<0.0001    *** P<0.000001\n\n\n');
fclose(fid) ;

%% 
x2 = yall(end,:);
se2 = seall(end,:);

z = (x1 - x2)./sqrt(se1.^2 + se2.^2);
p = normcdf(-z);

min(p);
max(p);

end
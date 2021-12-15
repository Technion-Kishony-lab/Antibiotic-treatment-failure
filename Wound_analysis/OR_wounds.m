function [] = OR_wounds(wound_cases, params)
% function takes the wound_cases structure and params and performs
% regresion for all susceptibitlyi-matched treated cases which gained resistance
% and remained sensitive. Fuction writes table of regressoin coefficients
% to the Tables directory.

any_previous_measurement = wound_cases.any_SRmeasurement;
previous_resistance = wound_cases.num_of_previous_R;
treatment_failures = wound_cases.treatfailure;
next_resistance_measurement = wound_cases.next_res;

%Adjusted for following demographics:
X_age = wound_cases.Demog.Age; X_age(:,6) = []; % reference age
zeroto39 = sum(X_age(:,1:4),2); %binned to larger bin for lower age groups since less data
eightyplus = sum(X_age(:,end-1:end),2);
X_demog = [zeroto39 X_age(:,5:end-2) eightyplus, wound_cases.Demog.Gender, wound_cases.Demog.Preg];
features = {'Age 0-39', 'Age 40-49', 'Age 60-69', 'Age 70-79', 'Age 80+', 'Gender', 'Pregnancy','Prev Resistance'};


min_1 = 5;
clear y_SSR p_SSR ciu_SSR cil_SSR se_SSR 
sen_res = params.resistant_group; % gain of resistance

X_all_to_use = zeros(size(any_previous_measurement(:,params.number_drugs)));
Y_all_to_use = zeros(size(any_previous_measurement(:,params.number_drugs)));
Y_all = zeros(size(any_previous_measurement(:,params.number_drugs)));
X_all = zeros(size(X_demog,1),size(X_demog,2)+1);


for drug = 1:params.number_drugs
SMP_to_use =  find(any_previous_measurement(:,drug)>0 & (wound_cases.SMP_Res(:,drug)== 1 | wound_cases.SMP_Res(:,drug)== 2) & next_resistance_measurement(:, drug) ~=0  );    
any_previous_R = previous_resistance(SMP_to_use ,drug)>0 ;
treatfailure = treatment_failures(SMP_to_use) ==1 & ismember(next_resistance_measurement(SMP_to_use, drug),sen_res) ;
drug_prch_temp = wound_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), any_previous_R(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
X_all(SMP_to_use(drug_prch_temp==1,:),:) = X;
X_all_to_use(SMP_to_use(drug_prch_temp==1,:)) = 1;
Y_all(SMP_to_use(drug_prch_temp==1)) =Y;
Y_all_to_use(SMP_to_use(drug_prch_temp==1,:)) = 1;
end

X_all(X_all_to_use == 0,:) = [];
Y_all(X_all_to_use == 0) = [];
c = my_fit(X_all, Y_all, min_1);
y_SSR(1) = exp(c.coef(end));
p_SSR(1) = c.p(end);
se_SSR(1) = c.se(end);
cil_SSR(1) = exp(c.coef(end) -  c.se(end));
ciu_SSR(1) = exp(c.coef(end) +  c.se(end));

%% print to table
features = pad(features);
yall = (c.coef(2:end));
pall = c.p(2:end);
seall = c.se(2:end);
CIL = (c.coef(2:end) -  c.se(2:end));
CIU = (c.coef(2:end) +  c.se(2:end));
pstrs = cell(size(pall)) ;
for i = 1:size(pall,1)
    for j = 1:size(pall,2)
        pstrs{i,j} = '' ;
        if pall(i,j)<0.01  , pstrs{i,j} = [pstrs{i,j}, '*']; end
        if pall(i,j)<0.001 , pstrs{i,j} = [pstrs{i,j}, '*']; end
        if pall(i,j)<0.0001, pstrs{i,j} = [pstrs{i,j}, '*']; end
    end
end
fid = fopen('Tables/Table_AdgustedOR_SR_wounds.txt', 'w');

    
fprintf(fid,'%-12s','');
fprintf(fid,'\n') ;
for ii = 1:size(yall,1)
    fprintf(fid,[features{ii} '  ']);
    coefs_tmp = [yall(ii) CIL(ii) CIU(ii) pstrs(ii)];
    fprintf(fid,'%5.2f [%5.2f,%5.2f]%-3s  ',coefs_tmp{:}) ;
    fprintf(fid,'%-7s','');
fprintf(fid,'\n') ;
end
fprintf(fid,'\n') ;
fprintf(fid,'* P<0.01    ** P<0.0001    *** P<0.000001\n\n\n');
fclose(fid) ;
%%
sen_res = params.sensitive_group; % remain sensitive

X_all_to_use = zeros(size(any_previous_measurement(:,drug)));
Y_all_to_use = zeros(size(any_previous_measurement(:,drug)));
Y_all = zeros(size(any_previous_measurement(:,drug)));
X_all = zeros(size(X_demog,1),size(X_demog,2)+1);

seall = zeros(size(X_demog,2) +2, params.number_drugs);
yall = zeros(size(X_demog,2) +2, params.number_drugs);
pall = zeros(size(X_demog,2) +2, params.number_drugs);

for drug = 1:params.number_drugs
SMP_to_use =  find(any_previous_measurement(:,drug)>0 & (wound_cases.SMP_Res(:,drug)== 1 | wound_cases.SMP_Res(:,drug)== 2) & next_resistance_measurement(:, drug) ~=0 );
any_previous_R = previous_resistance(SMP_to_use ,drug)>0 ;
treatfailure = treatment_failures(SMP_to_use) ==1 & ismember(next_resistance_measurement(SMP_to_use, drug),sen_res) ;
drug_prch_temp = wound_cases.PCR_sameday(SMP_to_use ,drug);
X_demog_to_use = X_demog(SMP_to_use,:);
X = [X_demog_to_use(drug_prch_temp==1,:), any_previous_R(drug_prch_temp==1)];
Y = treatfailure(drug_prch_temp==1);
X_all(SMP_to_use(drug_prch_temp==1,:),:) = X;
X_all_to_use(SMP_to_use(drug_prch_temp==1,:)) = 1;
Y_all(SMP_to_use(drug_prch_temp==1)) =Y;
Y_all_to_use(SMP_to_use(drug_prch_temp==1,:)) = 1;
end
X_all(X_all_to_use == 0,:) = [];
Y_all(X_all_to_use == 0) = [];
c = my_fit(X_all, Y_all, min_1);
y_SSR(2) = exp(c.coef(end));
p_SSR(2) = c.p(end);
se_SSR(2) = c.se(end);
cil_SSR(2) = exp(c.coef(end) -  c.se(end));
ciu_SSR(2) = exp(c.coef(end) +  c.se(end));

%% plot figure
b1 = bar(y_SSR, 'FaceColor','flat',  'BarWidth', 0.8);
hold on

b1(1).CData(1,:) = ones(1,1)*params.SR_color;
b1.CData(2,:) = ones(1,1)*params.SS_color;
%plot(b1.XData,cil_SSR,'kx');
%plot(b1.XData,ciu_SSR,'kx');
errneg = y_SSR - cil_SSR;
errpos =  ciu_SSR -y_SSR;

errorbar(b1.XData,y_SSR,errneg,errpos,'k', 'LineStyle','none');

hold on
for ii = 1:length(y_SSR)
    if p_SSR(ii)< 0.05 && p_SSR(ii)> 0.005
    plot(b1.XData(ii), y_SSR(ii)+1, '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.005 && p_SSR(ii)> 0.0005
    plot([b1.XData(ii)-0.1 b1.XData(ii)+0.1] , [y_SSR(ii)+1 y_SSR(ii)+1], '*k', 'MarkerSize',4)
    elseif p_SSR(ii) <= 0.0005
        plot([b1.XData(ii)-0.2  b1.XData(ii) b1.XData(ii)+0.2] , [ y_SSR(ii)+1  y_SSR(ii)+1 y_SSR(ii)+1], '*k', 'MarkerSize',4)
    end
end
% c.se(c.p> 0.05) = 0;
% hold on
% errorbar(c.coef(2:end),c.se(2:end),'.')
xticks((1:2:params.number_drugs*2)+0.5);
set(gca,'xticklabel','Wounds');
xtickangle(45);

ylabel({'Adjusted odds ratio of S-R failure for patients with';...
    'any previous resistance to the drug compared to sensitive'});

text(0,6.5,'{\it * P<0.05}');
text(0,6,'{\it** P<0.005}');
text(0,5.5,'{\it*** P<0.0005}');

set(gca, 'YScale', 'log');
%ylim([0.2 10])

%% print to table
features = pad(features);

yall = c.coef(2:end);
pall = c.p(2:end);
seall = c.se(2:end);
CIL = (c.coef(2:end) -  c.se(2:end));
CIU = (c.coef(2:end) +  c.se(2:end));

% yall = y_SSR(2:end,:);
% pall = p_SSR(2:end,:);
% seall= se_SSR(2:end,:);
% CIL = cil_SSR;
% CIU = ciu_SSR;

pstrs = cell(size(pall)) ;
for i = 1:size(pall,1)
    for j = 1:size(pall,2)
        pstrs{i,j} = '' ;
        if pall(i,j)<0.01  , pstrs{i,j} = [pstrs{i,j}, '*']; end
        if pall(i,j)<0.001 , pstrs{i,j} = [pstrs{i,j}, '*']; end
        if pall(i,j)<0.0001, pstrs{i,j} = [pstrs{i,j}, '*']; end
    end
end

fid = fopen('Tables/Table_AdgustedOR_SS_wounds.txt', 'w');

    
fprintf(fid,'%-12s','');
fprintf(fid,'\n') ;
for ii = 1:size(yall,1)
    fprintf(fid,[features{ii} '  ']);
    coefs_tmp = [yall(ii) CIL(ii) CIU(ii) pstrs(ii)];
    fprintf(fid,'%5.2f [%5.2f,%5.2f]%-3s  ',coefs_tmp{:}) ;
    fprintf(fid,'%-7s','');
fprintf(fid,'\n') ;
end
fprintf(fid,'\n') ;
fprintf(fid,'* P<0.01    ** P<0.0001    *** P<0.000001\n\n\n');
fclose(fid) ;

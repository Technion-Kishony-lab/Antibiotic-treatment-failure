for drug= 1:5
 relavant_suceptib_measured(:,drug) = wound_cases.PCR_sameday(:,drug) == 1  & wound_cases.SMP_Res(:, drug) ~= 0;
 relavant_suceptib_notmeasured(:,drug) = wound_cases.PCR_sameday(:,drug) == 1  & wound_cases.SMP_Res(:, drug) == 0;
 suceptib_matched(:,drug) =  wound_cases.PCR_sameday(:,drug) == 1  & ismember(wound_cases.SMP_Res(:, drug), [1 2]); 
 suceptib_mismatched(:,drug) =  wound_cases.PCR_sameday(:,drug) == 1  & ismember(wound_cases.SMP_Res(:, drug), 3); 
 treat_fail(:,drug) =  wound_cases.PCR_sameday(:,drug) == 1  & wound_cases.treatfailure; 
end

relavant_suceptib_measured_any = (sum(relavant_suceptib_measured,2));

nnz(relavant_suceptib_measured)

nnz(sum(suceptib_matched,2))
Run the main script UTI_analysis_and_figures.m to perform the analysis and generate the figures. 
The code requires the UTI_cases structure, which contains the following data:    

	  SMP_Res_drug_names: {8×1 cell}		% cell array of names of each of 8 antibiotic susceptibilties measured.
                    RandomID: [182571×1 double]		% array of unique patient identifier for each case
                     SMP_Res: [182571×8 double]		% array of current infection susceptibiltiy to each antibiotic. S = 1, I = 2, R = 3, 0 = unmeasured. If multiple isolates, maximal resistance is used.
                     hasdiag: [182571×1 logical]	% logical array is case had associated UTI diagnosis.
           num_of_previous_R: [182571×8 double] 	% array of number of previous infections patient has had which were resistant to each antibiotic.
           num_of_previous_S: [182571×8 double]		% array of number of previous infections patient has had which were sensitive to each antibiotic. 
           any_SRmeasurement: [182571×8 double]		% array of number of previous infections patient has had with any susceptibiltiy measured to each antibiotic 
        num_of_previous_SMPs: [182571×1 double]		% array of number of previous infections patient has had. 
                SamplingDate: [182571×1 double]		% array of date of current sample
                 changed_bac: [182571×1 double]		% logical array did future recurrent infection change species
                     new_bac: [182571×1 double] 	% array of species code for future recurrent infection.
                         bug: [182571×1 double]		% array of species code for current infection.
                     bug_all: {182571×1 cell}		% cell array of species code for all isolates in current infection.
                     new_bug: {182571×1 cell}		% cell array of species code for all isolates in future recurrent infection. NaN if no recurrence.
                        nBac: [182571×1 double]		% cell array of species code for all isolates in current infection.
                    next_res: [182571×8 double] 	% array of future recurrent infection maximal resistance to each antibiotic. NaN if no recurrence.
                         RES: {182571×8 cell}		% suceptibilty of each isolate in current infection.
                     new_RES: {182571×8 cell}		% suceptibilty of each isolate in future recurrent infection. NaN if no recurrence.
                   date_diff: [182571×1 double]		%number of days between current and future recurrent infection. NaN if no recurrence.
                treatfailure: [182571×1 logical]	% logical vector if any of the 5 antibiotics were prescribed 
    num_of_previous_R_byyear: [182571×8×4 double]	% array of number of previous infections patient has had which were resistant to each of 8 antibiotic for each 4 years prior.
    num_of_previous_S_byyear: [182571×8×4 double]	% array of number of previous infections patient has had which were sensitive to each of 8 antibiotic for each 4 years prior.
    any_SRmeasurement_byyear: [182571×8×4 double]	% array of number of previous infections patient has had with any susceptibiltiy measured to each of 8 antibiotic for each 4 years prior.
             PCR_sameday_any: [182571×1 logical]	% logical array if any of the 8 antibiotics were prescribed
                 PCR_sameday: [182571×9 logical]        % logical array of which antibitoic was prescribed. Column 9 is untreated.
           PCR_sameday_names: {1×9 cell}		% names of each prescribed antibitoics. 9 is untreated
                       Demog: [182571×8 table]   	% table of demographics of each case, including age bracket, gender, pregnancy, previous catheter use 
                        agem: [10×1 double] 	 	% midpoint of each age bracket in demograhics. 
              num_prev_purch: [182571×8 double]		% array of number of times the patient was previosuly prescribed each antibiotic.
                        Bugs: [609×2 table] 	 	% table of species names and codes.
                        
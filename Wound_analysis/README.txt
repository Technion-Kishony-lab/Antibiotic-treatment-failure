Run the main script wound_analysis_and_figures.m to perform the analysis and generate the figures. 
The code requires the wound_cases structure, which contains the following data:

      SMP_Res_drug_names: {5×1 cell} 		 % cell array of names of each antibiotic susceptibiltiy measured.
                RandomID: [7365×1 double]	 % array of unique patient identifier for each case
                 SMP_Res: [7365×5 double]	 % array of current infection susceptibiltiy to each antibiotic. S = 1, I = 2, R = 3, 0 = unmeasured. If multiple isolates, maximal resistance is used.
       num_of_previous_R: [7365×5 double] 	 % array of number of previous infections patient has had which were resistant to each antibiotic. 
       num_of_previous_S: [7365×5 double]	 % array of number of previous infections patient has had which were sensitive to each antibiotic. 
       any_SRmeasurement: [7365×5 double]	 % array of number of previous infections patient has had with any susceptibiltiy measured to each antibiotic 
    num_of_previous_SMPs: [7365×1 double]	 % array of number of previous infections patient has had. 
            SamplingDate: [7365×1 double]	 % array of date of current sample
                 bug_all: {7365×1 cell} 	 % cell array of species code for all isolates in current infection.
                 new_bug: {7365×1 cell} 	 % cell array of species code for all isolates in future recurrent infection. NaN if no recurrence.
                    nBac: [7365×1 double] 	 % number of isolated in current infection.
                     RES: {7365×5 cell} 	 % suceptibilty of each isolate in current infection.
                 new_RES: {7365×5 cell} 	 % suceptibilty of each isolate in future recurrent infection. NaN if no recurrence.
                next_res: [7365×5 double]	 % array of future recurrent infection maximal resistance to each antibiotic. NaN if no recurrence.
               date_diff: [7365×1 double] 	 % number of days between current and future recurrent infection. NaN if no recurrence.
            treatfailure: [7365×1 logical] 	 % logical array id case resulted in early recurrence within 28 days.
         PCR_sameday_any: [7365×1 logical] 	 % logical vector if any of the 5 antibiotics were prescribed 
             PCR_sameday: [7365×5 logical] 	 % logical array of which antibitoic was prescribed. 
       PCR_sameday_names: {1×5 cell} 	     % names of each prescribed antibitoics.
                   Demog: [7365×4 table] 	 % table of demographics of each wound case, including age bracket, gender, pregnancy. 
                    agem: [10×1 double] 	 % midpoint of each age bracket in demograhics. 
                    Bugs: [609×2 table] 	 % table of species names and codes.
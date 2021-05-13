README: 16S Analysis for Matsuda et al. 2021

In this folder you will find folders for each analysis run on these data:

For reference - SPECIESNAME refers to the coded name of the four species:
1. mcap = Montipora capitata
2. pcomp = Porites compressa
3. pvar = Pavona varians
4. pacu = Pocillopora acuta

ancom: analysis of compositions of microbiomes
	procedure: all R code for running ANCOMs
		> ancom_ace21.R = code for running ancoms
		> ANCOM_v2.1.R = code for ancom function
	output:
		by_phylotype: ANCOM contrasts by dominant algal phylotype
		by_treatment: ANCOM contrasts for each time point by treatment
			> ancom_SPECIESNAME_res.csv files contain all statistics (CLR mean difference & W stat) for each species & contrast
			> ancom_SPECIESNAME_figdata.csv files contain data for plotting (x = CLR mean diff, y = W stat)
			> SPECIESNAME_ancom.csv files have contain only taxa that are above cut-off of 60% (reject null hypothesis 60% of the time)
			used for bubble plots
			> SPECIESNAME_ancom.pdf volcano plots for each contrast (faceted T1 and TF for each species)
			
bac_figures: Code for combined microbiome metric bubble plots
	procedure: figures_ace21.R = code for creating figures
	output: 
		> SPECIESNAME_top_taxa.csv files contain combined abundances of "top" taxa: top most abundant, most shared and most differential
		> SPECIESNAME_top_bubble.pdf files contain figures (also saved as png if pdf doesn't open) 

core_microbiome: Shared taxa 
	procedure: core_ace21.R
	output: SPECIESNAME_core.csv files contain all core taxa for each species based on varying percentages
		> mcap = 60%, pcomp = 50%, pvar = 45% and pacu = 57% (see methods)
		
indicator_analysis: multi-level pattern analysis
	procedure: indicator_analysis_ace21.R 
	output: 
		> species_inidicval_with_tax.csv contains results from indicator analysis in table form, with concatenated taxa names
		> species_indicval.csv file has all data without taxonomy

top_abundance: Top most abundant taxa
	procedure: abundance_ace21.R
	output: SPECIESNAME_abund.csv files contain the top 10 most abundant taxa for each species (all samples)










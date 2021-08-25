# Single_SNP_Matching_Transplant
Code for data analysis and statistical simulations for single SNP matching methods for paired D/R transplant data

Results are available in published manuscript at https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.22349

Citation:
Arthur, VL, Guan, W, Loza, B-l, Keating, B, Chen, J. Joint testing of donor and recipient genetic matching scores and recipient genotype has robust power for finding genes associated with transplant outcomes. Genetic Epidemiology. 2020; 44: 893– 907. https://doi.org/10.1002/gepi.22349

Directories:

-Paper1_CodeToSendToCollabs: Code used to conduct real data analysis on time-to-event data for liver transplant. Code was sent to collaborators to conduct the testing.

-Gene_Files: Files with SNP info from the genes used in simulation studies. All data is from publically available HapMap3 datasets.

-VCFheaderfiles: Headers for creating VCF files when running simulations.

-Related_Pairs_Files: File to pair simulated haplotypes into genotype data

-Create_Tables_And_Figures: Code used to create all tables and figures in the manuscript

-All_Scripts_For_Paper1_Simulation_Pipeline: All scripts used in the simulation pipeline for the manuscript

----------------------------------------------------------
Simulation Pipeline/Order of Scripts

*For Type I Error Analysis and Power Analysis
1) SingleGenePipeline_NormBeagle.sh
   Required inputs: 
		- Chromosome Number for Gene where Haplotypes are coming from
		- Sample Size (Number of D/R pairs)
		- Number of Simulations to Run
		- Random Number (seed)
		- Name of Gene where Haplotypes are from
		- Which Simulation to Start at (allows for running clusters of sims)

	Dependencies:
	- PullHaplotypesForSimulations.sh
		- Requires R v 3.6.1
		- Haplotype Data Files (plink format)
	- AddVCFHeaders.sh
		- VCF Header File for sample size
	- ConvertVCFToPlink.sh
		- Requires Plink v 1.9
	- ConvertPlinkToBeagle.sh
		- Requires Plink v 1.9
	- RunBeagleIBD.sh
		- Requires Java and Beagle v 3.3.2
		- File specifying related pairs (D/R pairs)
	- CalcIBSAndIncompScoresVariableOnly.sh
		- Requires R v 3.5.3
			Installed Packages: ARTP2, dplyr
	- ReformatBeagleIBDScores.sh
		- Requires R v 3.5.3
	- CalcPValues_TypeIErr.sh 
		- Requires R v 3.5.3
	- CalcPValues_PowerAnalysis.sh
		- Requires R v 3.5.3
	
	Outputs:
		- Files with p-values for both type I error and power analyses
		
2) Calc_TypeIErr.sh
   Required Inputs:
		- Chromosome number of Gene where Haplotypes were sampled from
		- Samples Size (Number of D/R pairs)
		- Number of Simulations Run in Step 1
		- Name of Gene where Haplotypes were sampled from
	
	Dependencies:
	- Requires R v 3.6.1
	- CalcTypeIError.R
	
3) Calc_Power.sh
	Required Inputs:
		- Chromosome number of Gene where Haplotypes were sampled from
		- Samples Size (Number of D/R pairs)
		- Number of Simulations Run in Step 1
		- Name of Gene where Haplotypes were sampled from
		
	Dependencies:
	- Requires R v 3.6.1
	- CalcPower.R

4) TIE_plots_Haplotypes.R (Run Locally to plot Type I Error)
	Dependencies:
	- Requires R v 3.6.1
		Installed Packages: RColorBrewer
		
5) Power_BarCharts_TNFAIP8L3.R (Run Locally to plot Marginal Power)
	Dependencies:
	- Requires R v 3.6.1
		Installed Packages: ggplot2, RColorBrewer
		
*For joint power analyses
1) Need to run SingleGenePipeline_NormBeagle.sh first to obtain calculated score .csv files

2) JointPowerAnalyses.sh 
	Required inputs: 
		- Chromosome Number for Gene where Haplotypes are from
		- Sample Size (Number of D/R pairs)
		- Number of Simulations to Run
		- Random Number (seed)
		- Name of Gene where Haplotypes are from
		- Which Simulation to Start at (allows for running clusters of sims)
	
	Dependencies:
		- Run_JointlyTestingPowerScoreRGeno.sh
			- Requires R v 3.5.3
			- JointlyTestingPowerScoreRGeno.R
				Required Packages: epicalc

3) Run_CalcJointPower.sh
	Required Inputs: 
		- Chromosome Number for Gene where Haplotypes are from
		- Sample Size (Number of D/R pairs)
		- Number of Simulations to Run
		- Name of Gene where Haplotypes are from
		
	Dependencies:
	- Requires R 3.5.3
	- CalcJointPower.R
		
4) JointPowerAnalyses_RGenoTrue.sh
	Required inputs: 
		- Chromosome Number for Gene where Haplotypes are from
		- Sample Size (Number of D/R pairs)
		- Number of Simulations to Run
		- Random Number (seed)
		- Name of Gene where Haplotypes are from
		- Which Simulation to Start at (allows for running clusters of sims)
		
	Dependencies:
		- Run_JointlyTestingPowerRGenoScore.sh
			- Requires R v 3.5.3
			- JointlyTestingPowerRgenoScore.R
				Required Packages: epicalc
				
5) Run_CalcJointPower_RGenoTrue.sh
Required Inputs: 
		- Chromosome Number for Gene where Haplotypes are from
		- Sample Size (Number of D/R pairs)
		- Number of Simulations to Run
		- Name of Gene where Haplotypes are from
		
	Dependencies:
	- Requires R 3.5.3
	- CalcJointPower_RGenoTrue.R
	
6) Plotting_JointPowerTests_TNFAIP8L3.R (Run locally to plot all joint power test results)
	-Requires R v 3.6.1
		Required Packages: RColorBrewer, ggplot2

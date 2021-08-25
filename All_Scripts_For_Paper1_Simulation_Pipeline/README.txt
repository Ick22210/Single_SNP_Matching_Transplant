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
				

				
		
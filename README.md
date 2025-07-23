# Long-term effects of sweeteners and sweetness enhancers consumption after weight loss on the gut microbiota composition in individuals with overweight or obesity: The SWEET study

This is a complete microbiota and primary outcomes analysis pipeline. To reproduce all results and figures, first use the bash script to generate the ASV table and associated objects with the qiime2. Then, use qiime2 artifacts to make a phyloseq object and run the script in order. 

## Abstract 

**Background and objectives:** The long-term effects of sweeteners and sweetness enhancers (S&SEs) on the human gut microbiome remain elusive, with prior studies showing inconsistent results. In the multi-centre SWEET project, we investigated the effect of S&SEs as a replacement of sugar on gut microbial composition and safety outcomes in adults with overweight or obesity.
<br />
<br />
**Methods:** In this randomized controlled trial, 341 participants (age 18-65 yrs and BMI≥25 kg/m2) from different European countries (The Netherlands, Greece, Spain, and Denmark) underwent a 2-month weight loss (WL) phase, followed by a 10-month weight maintenance (WM) phase during which they adhered to a healthy ad libitum diet (<10 energy % added sugar) either with or without S&SEs products (S&SEs vs. sugar group). Fecal samples were analyzed using 16S rRNA sequencing in a subgroup of 137 completers at baseline (month M0), after WL (M2), M6, and M12 to determine gut microbial composition. Liver fat content was assessed in a subgroup of 29 participants using proton-magnetic resonance spectroscopy at M0, M2, and M12. Adverse events and changes in medication use were monitored in all participants throughout the study.
<br />
<br />
**Results:** The S&SEs group demonstrated lower weight regain (change M12-M2: 3.4±0.7 vs. 5.6±0.8 kg, P=0.011) and lower energy intake (P=0.044) compared to the sugar group. Distinct shifts in microbial composition were observed between groups, with a higher abundance of taxa related to short-chain fatty acid (SCFA) and methane production in the S&SEs compared to the sugar group (q≤0.05). The S&SEs group reported more gastrointestinal symptoms, including abdominal pain, loose stools, and excess gas than the sugar group. Changes in liver fat and concomitant medication use were not different between groups. Baseline microbiota composition and/or composition after initial weight loss could classify response to intervention with respect to body weight regain and change in HbA1c with reasonable accuracy.
<br />
<br />
**Discussion:** The S&SEs group altered gut microbial composition towards a higher abundance of SCFA and methane-producing bacterial taxa compared to the sugar group during a 10-month WM phase, which was accompanied by an improved body weight control, no changes in cardiometabolic health and more gastrointestinal symptoms.  Microbial composition could classify response and non-response to S&SE with reasonable accuracy, suggesting that effects may be to some extent person specific. 

# Results reproduction

This repository contains all data and scripts to reproduce all figures and tables (limited formatting options).  To reproduce results, clone the repository and run "run_all.R". Alternatively, scripts can be run separately; however, firstly, the "prm.R" and then "0_data_prep.R" should be executed.  
<br />
The information about package versions can be found in the "sessionInfo.txt" file. In addition, the files for the automatic setting of the environment with "renv" package are provided. 
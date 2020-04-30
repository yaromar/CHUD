## Maryam Zekavat

R

library(data.table)
pheno7 = fread("/path/to/phenoFile.gz",header=T)
pheno7 = data.frame(pheno7)

varcounts = data.frame()
for (i in 1:10){
tmp = fread(paste('/path/to/CHUD/VarCounts/varcounts', i, '.txt', sep=""), header=TRUE, sep="\t")
tmp = data.frame(tmp)
varcounts = rbind(varcounts,tmp)
}
varcounts$sample_id = unlist(lapply( strsplit(varcounts$sample_id,"_"), "[", 1))

pheno8 = merge(pheno7, varcounts, by=c(1))
pheno8$Large_CHIP = factor(pheno8$Large_CHIP , levels=c("NO_CHIP", "LARGE_CHIP","SMALL_CHIP"))
vars=colnames(varcounts)[2:18]
pheno_list2 = c("AML", "MPN")
vars = c(vars, c('hasCHIP', 'Large_CHIP'))
summaryDF = data.frame()
tmp=pheno8

## Running logistic regression:
for (i in 1:length(pheno_list2)){
for (j in 1:length(vars)){

print(pheno_list2[i])
indx = which(colnames(tmp) == vars[j])[1]
pheno_col = which(colnames(tmp) ==pheno_list2[i])


df = as.data.frame(t(as.data.frame(summary( glm(tmp[,pheno_col] ~ tmp[,indx]+age +age2 + Sex_numeric+ PC1+PC2+PC3+PC4+PC5, data = tmp, family='binomial'))$coeff[2,])))
df$x = vars[j]
df$y = pheno_list2[i]
df$N_ttl = length(which(!is.na(tmp[,indx]) & !is.na(tmp[,pheno_col])& !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$ever_smoked)))
df$N_Cases = length(which(tmp[,pheno_col] == 1 & !is.na(tmp[,indx]) & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$ever_smoked)))
df$N_Controls = length(which(tmp[,pheno_col] == 0  & !is.na(tmp[,indx]) & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$ever_smoked)))
df$N_Cases_withVar = length(which(tmp[,pheno_col] == 1 & (tmp[,indx] >0 | tmp[,indx] == "LARGE_CHIP") & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$ever_smoked)))
df$N_Controls_withVar = length(which(tmp[,pheno_col] == 0  &(tmp[,indx] >0 | tmp[,indx] == "LARGE_CHIP")  & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$ever_smoked)))
df$adjustment = "adjusted"

summaryDF= rbind(summaryDF, df)

df = as.data.frame(t(as.data.frame(summary( glm(tmp[,pheno_col] ~ tmp[,indx], data = tmp, family='binomial'))$coeff[2,])))
df$x = vars[j]
df$y = pheno_list2[i]
df$N_ttl = length(which(!is.na(tmp[,indx]) & !is.na(tmp[,indx]) ))
df$N_Cases = length(which(tmp[,pheno_col] == 1 & !is.na(tmp[,indx]) ))
df$N_Controls = length(which(tmp[,pheno_col]== 0  & !is.na(tmp[,indx]) ))
df$N_Cases_withVar = length(which(tmp[,pheno_col] > 0 & (tmp[,indx] >0 | tmp[,indx] == "LARGE_CHIP")))
df$N_Controls_withVar = length(which(tmp[,pheno_col] == 0 & (tmp[,indx] >0 | tmp[,indx] == "LARGE_CHIP") ))
df$adjustment = "not adjusted"
summaryDF = rbind(summaryDF, df)

}
}
write.table(summaryDF,"/medpop/esp2/mzekavat/CHIP/CHUD/VarCounts.PheWAS.logreg.txt",
col.names = T, row.names = F, quote = F, sep = "\t")

#Separately running associations adjusted for CHIP:
summary( glm(AML ~overlap_._LeukemiaGene+hasCHIP + age +age2 + Sex_numeric+ PC1+PC2+PC3+PC4+PC5, data = pheno8, family='binomial'))
Coefficients:
                         Estimate Std. Error z value Pr(>|z|)    
(Intercept)            -15.045889  10.035398  -1.499    0.134    
overlap_._LeukemiaGene   0.514402   0.123809   4.155 3.26e-05 ***
hasCHIP                  0.752929   0.462734   1.627    0.104    
age                      0.229124   0.347006   0.660    0.509    
age2                    -0.001441   0.002992  -0.482    0.630    
Sex_numeric              0.457614   0.298849   1.531    0.126    
PC1                      0.097978   0.094445   1.037    0.300    
PC2                      0.142359   0.098202   1.450    0.147    
PC3                     -0.025322   0.095088  -0.266    0.790    
PC4                     -0.083292   0.065686  -1.268    0.205    
PC5                      0.011948   0.032049   0.373    0.709    

summary( glm(AML ~overlap_._LeukemiaGene+ age +age2 + Sex_numeric+ PC1+PC2+PC3+PC4+PC5, 
		data = pheno8[-which(pheno8$hasCHIP==1),], family='binomial'))
Coefficients:
                        Estimate Std. Error z value Pr(>|z|)    
(Intercept)            -13.28749   10.13036  -1.312 0.189638    
overlap_._LeukemiaGene   0.48392    0.13771   3.514 0.000441 ***
age                      0.18030    0.35264   0.511 0.609157    
age2                    -0.00105    0.00306  -0.343 0.731411    
Sex_numeric              0.53112    0.32168   1.651 0.098724 .  
PC1                      0.08319    0.10075   0.826 0.409002    
PC2                      0.02771    0.10495   0.264 0.791758    
PC3                     -0.04441    0.10182  -0.436 0.662707    
PC4                     -0.09971    0.06925  -1.440 0.149904    
PC5                      0.01961    0.03362   0.583 0.559810    


summary( glm(AML ~overlap_._LeukemiaGene_._binom_._VAF+ age +age2 + Sex_numeric+ PC1+PC2+PC3+PC4+PC5, 
		data = pheno8[-which(pheno8$hasCHIP==1 ),], family='binomial'))
Coefficients:
                                       Estimate Std. Error z value Pr(>|z|)    
(Intercept)                          -13.285536  10.160303  -1.308    0.191    
overlap_._LeukemiaGene_._binom_._VAF   1.203503   0.300619   4.003 6.24e-05 ***
age                                    0.184902   0.353788   0.523    0.601    
age2                                  -0.001103   0.003069  -0.359    0.719    
Sex_numeric                            0.518037   0.321830   1.610    0.107    
PC1                                    0.083357   0.101152   0.824    0.410    
PC2                                    0.029513   0.104949   0.281    0.779    
PC3                                   -0.046850   0.101494  -0.462    0.644    
PC4                                   -0.098732   0.069367  -1.423    0.155    
PC5                                    0.019228   0.033580   0.573    0.567    

summary( glm(AML ~overlap_._LeukemiaGene_._binom_._VAF+ age +age2 + Sex_numeric+ PC1+PC2+PC3+PC4+PC5, 
		data = pheno8[-which(pheno8$hasCHIP==1 | pheno8$overlap_._TOPMed_CHIPVar ==1),], family='binomial'))

Coefficients:
                                       Estimate Std. Error z value Pr(>|z|)    
(Intercept)                          -16.751608  10.820213  -1.548 0.121580    
overlap_._LeukemiaGene_._binom_._VAF   1.184322   0.332021   3.567 0.000361 ***
age                                    0.315768   0.377881   0.836 0.403364    
age2                                  -0.002324   0.003288  -0.707 0.479756    
Sex_numeric                            0.367802   0.333497   1.103 0.270085    
PC1                                    0.077548   0.106472   0.728 0.466403    
PC2                                    0.015033   0.110194   0.136 0.891488    
PC3                                   -0.114759   0.106412  -1.078 0.280835    
PC4                                   -0.121273   0.072835  -1.665 0.095906 .  
PC5                                    0.030784   0.034785   0.885 0.376170    

summary( glm(AML ~overlap_._SOMATIC+ age +age2 + Sex_numeric+ PC1+PC2+PC3+PC4+PC5, 
		data = pheno8[-which(pheno8$hasCHIP==1),], family='binomial'))
Coefficients:
                    Estimate Std. Error z value Pr(>|z|)    
(Intercept)       -13.323213  10.158149  -1.312 0.189662    
overlap_._SOMATIC   0.041182   0.011714   3.516 0.000439 ***
age                 0.178351   0.353481   0.505 0.613870    
age2               -0.001014   0.003066  -0.331 0.740913    
Sex_numeric         0.525455   0.321653   1.634 0.102341    
PC1                 0.086220   0.100802   0.855 0.392363    
PC2                 0.034032   0.105403   0.323 0.746792    
PC3                -0.044821   0.101682  -0.441 0.659364    
PC4                -0.100907   0.069269  -1.457 0.145188    
PC5                 0.021489   0.033605   0.639 0.522520    

summary( glm(AML ~overlap_._SOMATIC_._binom_._VAF+ age +age2 + Sex_numeric+ PC1+PC2+PC3+PC4+PC5, 
		data = pheno8[-which(pheno8$hasCHIP==1),], family='binomial'))
Coefficients:
                                  Estimate Std. Error z value Pr(>|z|)  
(Intercept)                     -13.295810  10.159758  -1.309   0.1906  
overlap_._SOMATIC_._binom_._VAF   0.247206   0.122516   2.018   0.0436 *
age                               0.178820   0.353411   0.506   0.6129  
age2                             -0.001024   0.003066  -0.334   0.7384  
Sex_numeric                       0.524196   0.321505   1.630   0.1030  
PC1                               0.081175   0.101037   0.803   0.4217  
PC2                               0.032572   0.105475   0.309   0.7575  
PC3                              -0.047016   0.101712  -0.462   0.6439  
PC4                              -0.100225   0.069570  -1.441   0.1497  
PC5                               0.020821   0.033704   0.618   0.5367  

summary( glm(AML ~overlap_._SOMATIC+ age +age2 + Sex_numeric+ PC1+PC2+PC3+PC4+PC5, 
		data = pheno8[-which(pheno8$hasCHIP==1 | pheno8$overlap_._TOPMed_CHIPVar ==1),], family='binomial'))
Coefficients:
                    Estimate Std. Error z value Pr(>|z|)   
(Intercept)       -16.871901  10.818992  -1.559  0.11889   
overlap_._SOMATIC   0.040311   0.012952   3.112  0.00186 **
age                 0.313131   0.377648   0.829  0.40701   
age2               -0.002277   0.003286  -0.693  0.48840   
Sex_numeric         0.375167   0.333216   1.126  0.26021   
PC1                 0.081217   0.105892   0.767  0.44309   
PC2                 0.018469   0.110757   0.167  0.86756   
PC3                -0.114054   0.106909  -1.067  0.28605   
PC4                -0.122472   0.072550  -1.688  0.09139 . 
PC5                 0.032668   0.034749   0.940  0.34715   


## Running associations with age
vars = c("age")
pheno_list2=c(colnames(varcounts)[2:18], c('hasCHIP', 'Large_CHIP'))
summaryDF = data.frame()

for (i in 1:length(pheno_list2)){
for (j in 1:length(vars)){

print(pheno_list2[i])
indx = which(colnames(tmp) == vars[j])[1]
pheno_col = which(colnames(tmp) ==pheno_list2[i])

df = as.data.frame(t(as.data.frame(summary( glm(tmp[,pheno_col] ~ tmp[,indx], data = tmp))$coeff[2,])))
df$x = vars[j]
df$y = pheno_list2[i]
df$N_ttl = length(which(!is.na(tmp[,indx]) & !is.na(tmp[,indx]) ))
df$N_withVar = length(which(tmp[,pheno_col] == 1 &  (tmp[,indx] >0 | tmp[,indx] == "LARGE_CHIP")  ))
summaryDF = rbind(summaryDF, df)

}
}

write.table(summaryDF,"/medpop/esp2/mzekavat/CHIP/CHUD/VarCounts.age.linreg.txt",
col.names = T, row.names = F, quote = F, sep = "\t")

### Running CoxPh associations with CAD
library(survival)
pheno_list = c("Coronary_Artery_Disease_SOFT")
vars=c(colnames(varcounts)[2:18], c('hasCHIP', 'Large_CHIP'))

summaryDF=data.frame()
for (i in 1:length(pheno_list)){
	print(pheno_list[i])
for (j in 1:length(vars)){
	print(vars[j])
tmp = pheno8
indx = which(colnames(tmp) == vars[j])[1]

#removing prevalent cases
preval_col = which(colnames(tmp) == paste("Prev_",pheno_list[i],sep=""))
incid_col = which(colnames(tmp) == paste("Incd_",pheno_list[i],sep=""))

tmp = tmp[-which(tmp[,preval_col] == 1),]

if (length(which(tmp[,incid_col] == 1 & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$ever_smoked))) > 10){

#running survival anal
tmp$SurvObj <- Surv(tmp[,which(colnames(tmp) == paste(pheno_list[i], "_FollowUp", sep=""))],  tmp[,which(colnames(tmp) == paste("Incd_",pheno_list[i], sep=""))]== 1)

df = as.data.frame(t(as.data.frame(summary( coxph(tmp$SurvObj ~ tmp[,indx]+age +age2 + Sex_numeric+PC1+PC2+PC3+PC4+PC5, data = tmp))$coeff[1,])))
df$y = pheno_list[i]
df$x = vars[j]
df$N_Incd_cases = length(which(tmp[,incid_col] == 1 & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$ever_smoked)))
df$N_Controls = length(which(tmp[,incid_col] == 0 & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$ever_smoked)))
df$N_Incd_cases_withVar = length(which(tmp[,incid_col] ==1 & (tmp[,indx] >0 | tmp[,indx] == "LARGE_CHIP") & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$ever_smoked)))
df$N_Controls_withVar = length(which(tmp[,incid_col] == 0 & (tmp[,indx] >0 | tmp[,indx] == "LARGE_CHIP") & !is.na(tmp$age) & !is.na(tmp$PC1)  & !is.na(tmp$ever_smoked)))
summaryDF = rbind(summaryDF, df)
}
}
}
write.table(summaryDF,"/medpop/esp2/mzekavat/CHIP/CHUD/VarCounts.cad_coxph.txt",
col.names = T, row.names = F, quote = F, sep = "\t")


##### Making forest plots:
library(ggplot2) #plotting system
library(grid)
library(reshape2)
library(scales)
library(ggrepel)
library(plotly)
library(meta)
Forest = read.table("/medpop/esp2/mzekavat/CHIP/CHUD/data/Forest_logreg_AML_MPN_CHUD.txt" ,
                   header=T, as.is=T, stringsAsFactors=F, comment.char = '', sep="\t")
Forest$Pval = as.character(Forest$Pr...z..)
Forest$x = ifelse(Forest$x == "overlap_._TOPMed_CHIPVar_._VAF", 'Rare_deleterious_TOPMed_CHIPvar_LargeClone',
				ifelse(Forest$x == "overlap_._LeukemiaGene_._binom_._VAF", 'Rare_deleterious_LeukemiaGene_binom_LargeClone', 
					ifelse(Forest$x == "overlap_._TOPMed_CHIPVar", 'TOPMed_CHIPvar',
						ifelse(Forest$x == "overlap_._LeukemiaGene_._VAF", "Rare_deleterious_LeukemiaGene_LargeClone",
							ifelse(Forest$x == "overlap_._TOPMed_CHIPVar_._binom_._VAF", 'Rare_deleterious_TOPMed_CHIPvar_binom_LargeClone',
								ifelse(Forest$x == "overlap_._TOPMed_CHIPVar_._binom", 'Rare_deleterious_TOPMed_CHIPvar_binom',
									ifelse(Forest$x  == "overlap_._LeukemiaGene", "Rare_deleterious_LeukemiaGene",
										ifelse(Forest$x == "overlap_._LeukemiaGene_._binom", "Rare_deleterious_LeukemiaGene_binom",
											ifelse(Forest$x == "overlap_._binom","Rare_deleterious_binom",
											ifelse(Forest$x == "overlap", "Rare_deleterious",
												ifelse(Forest$x == "overlap_._SOMATIC", "Rare_deleterious_SOMATIC",
													ifelse(Forest$x == "overlap_._SOMATIC_._binom","Rare_deleterious_SOMATIC_binom",
														ifelse(Forest$x == "num_FilterMutect","All_QCed_Mutect2",
															ifelse(Forest$x == "overlap_._SOMATIC_._binom_._VAF",'Rare_deleterious_SOMATIC_binom_LargeClone',
																ifelse(Forest$x == "overlap_._binom_._VAF", "Rare_deleterious_binom_LargeClone", 
																	ifelse(Forest$x == "overlap_._VAF", "Rare_deleterious_LargeClone",
																		ifelse(Forest$x == "overlap_._SOMATIC_._VAF","Rare_deleterious_SOMATIC_LargeClone", Forest$x )))))))) )))))))))

summarymeta=metagen(Estimate,Std..Error,studlab=x,byvar = y,data=Forest,sm="OR")
pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/Forest_logreg_AML_MPN_CHUD.pdf",sep=""), width = 18, height= 10)
forest(summarymeta,leftcols =c("studlab"),leftlabs = c(""),just='left',colgap.left=unit(0.1,"cm"),xlim=c(1,12),at=c(1,3,6,9,12),
      smlab='Effect of Somatic Variant (N) \non Hematologic Cancer',colgap=unit(7, "mm"),rightcols=c("OR",'ci','Pval', 'N_Cases','N_Controls','N_Cases_withVar','N_Controls_withVar'),
      rightlabs=c("OR",'95% CI','P', "Cases (N)", "Controls (N)", "Cases with Var>0 (N)", "Controls with Var>0 (N)"),plotwidth=unit(6, "cm"),
      comb.random=F,print.Q=F,overall=F,comb.fixed=F,print.byvar=F,addspace=T,xlab = "OR")
dev.off()

Forest = read.table("/medpop/esp2/mzekavat/CHIP/CHUD/data/Forest_hr_CAD_CHUD.txt" ,
                   header=T, as.is=T, stringsAsFactors=F, comment.char = '', sep="\t")
Forest$Pval = as.character(Forest$P)
Forest$x = ifelse(Forest$x == "overlap_._TOPMed_CHIPVar_._VAF", 'Rare_deleterious_TOPMed_CHIPvar_LargeClone',
				ifelse(Forest$x == "overlap_._LeukemiaGene_._binom_._VAF", 'Rare_deleterious_LeukemiaGene_binom_LargeClone', 
					ifelse(Forest$x == "overlap_._TOPMed_CHIPVar", 'TOPMed_CHIPvar',
						ifelse(Forest$x == "overlap_._LeukemiaGene_._VAF", "Rare_deleterious_LeukemiaGene_LargeClone",
							ifelse(Forest$x == "overlap_._TOPMed_CHIPVar_._binom_._VAF", 'Rare_deleterious_TOPMed_CHIPvar_binom_LargeClone',
								ifelse(Forest$x == "overlap_._TOPMed_CHIPVar_._binom", 'Rare_deleterious_TOPMed_CHIPvar_binom',
									ifelse(Forest$x  == "overlap_._LeukemiaGene", "Rare_deleterious_LeukemiaGene",
										ifelse(Forest$x == "overlap_._LeukemiaGene_._binom", "Rare_deleterious_LeukemiaGene_binom",
											ifelse(Forest$x == "overlap_._binom","Rare_deleterious_binom",
											ifelse(Forest$x == "overlap", "Rare_deleterious",
												ifelse(Forest$x == "overlap_._SOMATIC", "Rare_deleterious_SOMATIC",
													ifelse(Forest$x == "overlap_._SOMATIC_._binom","Rare_deleterious_SOMATIC_binom",
														ifelse(Forest$x == "num_FilterMutect","All_QCed_Mutect2",
															ifelse(Forest$x == "overlap_._SOMATIC_._binom_._VAF",'Rare_deleterious_SOMATIC_binom_LargeClone',
																ifelse(Forest$x == "overlap_._binom_._VAF", "Rare_deleterious_binom_LargeClone", 
																	ifelse(Forest$x == "overlap_._VAF", "Rare_deleterious_LargeClone",
																		ifelse(Forest$x == "overlap_._SOMATIC_._VAF","Rare_deleterious_SOMATIC_LargeClone", Forest$x )))))))) )))))))))

summarymeta=metagen(coef,se.coef.,studlab=x,data=Forest,sm="HR")
pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/Forest_hr_CAD_CHUD.pdf",sep=""), width = 18, height= 6)
forest(summarymeta,leftcols =c("studlab"),leftlabs = c(""),just='left',colgap.left=unit(0.1,"cm"),xlim=c(.95,2.5),at=c(.95,1,1.1,1.5,2,2.5),
      smlab='Effect of Somatic Variant (N) \non CAD',colgap=unit(7, "mm"),rightcols=c("HR",'ci','Pval', 'N_Incd_cases','N_Controls','N_Incd_cases_withVar','N_Controls_withVar'),
      rightlabs=c("HR",'95% CI','P', "Incd Cases (N)", "Controls (N)", "Incd Cases with Var>0 (N)", "Controls with Var>0 (N)"),plotwidth=unit(6, "cm"),
      comb.random=F,print.Q=F,overall=F,comb.fixed=F,print.byvar=F,addspace=T,xlab = "HR")
dev.off()


#### Making figures
library(tidyr)

pheno9= pheno8[,c(1,3,4,278, 640,644,631,632, 633,720:735)]
pheno9$Large_CHIP = ifelse(pheno9$Large_CHIP == "LARGE_CHIP", 1, 0)
> colnames(pheno9)
 [1] "id"                                    
 [2] "Sex_numeric"                           
 [3] "age"                                   
 [4] "Coronary_Artery_Disease_SOFT"          
 [5] "AML"                                   
 [6] "MPN"                                   
 [7] "AF"                                    
 [8] "hasCHIP"                               
 [9] "Large_CHIP"                            
[10] "overlap"                               
[11] "overlap_._LeukemiaGene"                
[12] "overlap_._TOPMed_CHIPVar"              
[13] "overlap_._SOMATIC"                     
[14] "overlap_._binom"                       
[15] "overlap_._LeukemiaGene_._binom"        
[16] "overlap_._TOPMed_CHIPVar_._binom"      
[17] "overlap_._SOMATIC_._binom"             
[18] "overlap_._VAF"                         
[19] "overlap_._LeukemiaGene_._VAF"          
[20] "overlap_._TOPMed_CHIPVar_._VAF"        
[21] "overlap_._SOMATIC_._VAF"               
[22] "overlap_._binom_._VAF"                 
[23] "overlap_._LeukemiaGene_._binom_._VAF"  
[24] "overlap_._TOPMed_CHIPVar_._binom_._VAF"
[25] "overlap_._SOMATIC_._binom_._VAF"       

colnames(pheno9)[10:25] = c('Rare_deleterious', 'Rare_deleterious_LeukemiaGene', 'TOPMed_CHIPvar', 'Rare_deleterious_SOMATIC', 
							'Rare_deleterious_binom', 'Rare_deleterious_LeukemiaGene_binom', 'Rare_deleterious_TOPMed_CHIPvar_binom', 'Rare_deleterious_SOMATIC_binom',
							'Rare_deleterious_LargeClone','Rare_deleterious_LeukemiaGene_LargeClone', 'Rare_deleterious_TOPMed_CHIPvar_LargeClone', 'Rare_deleterious_SOMATIC_LargeClone',
							'Rare_deleterious_binom_LargeClone', 'Rare_deleterious_LeukemiaGene_binom_LargeClone', 'Rare_deleterious_TOPMed_CHIPvar_binom_LargeClone','Rare_deleterious_SOMATIC_binom_LargeClone')

pheno9 = pheno9[,c(1:7,10:25,8,9)]
counts_and_phenos_LONG = gather(pheno9, Filtration, VarCount, Rare_deleterious:Large_CHIP, factor_key=TRUE)



binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

#subset$Filtration = factor(subset$Filtration , levels=c("SOMATIC", "LeukemiaGene","hasCHIP"))

p <- ggplot(counts_and_phenos_LONG, aes(x = age, y = VarCount, group = Filtration, color=as.factor(Filtration)))
p <- p +geom_smooth()#+ylim(0,0.4)#method = lm, formula = y ~ splines::bs(x, 3))+ylim(0,0.4)
p <- p +facet_wrap(~Filtration, scales="free",ncol=4)+  xlab("Age") + ylab("Prevalence")
p <- p + theme(
	strip.background = element_blank(),
	strip.text.x = element_text(size = 14),
	strip.text.y = element_text(size = 14),
	axis.text.x = element_text(size=13, hjust=1, color = "black"),
	axis.text.y = element_text(size=13,hjust=1, color = "black"),
	axis.ticks =  element_line(colour = "black"), 
	axis.title.y= element_text(size=14),
	axis.title.x= element_text(size=14),
	legend.position="right", 
	legend.title = element_blank(),
	panel.background = element_blank(), 
	panel.border = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	panel.spacing = unit(1.0, "lines"), 
	plot.background = element_blank(), 
	plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
	axis.line.x = element_line(colour = "black"),
	axis.line.y = element_line(colour = "black")
)+theme(legend.position="none")

pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/LeukemiaVar.ageAssoc_spline.v2.pdf",sep=""), width = 18, height= 18)
print(p)
dev.off()

subset = counts_and_phenos_LONG[which(counts_and_phenos_LONG$Filtration %in% c('Rare_deleterious_binom_LargeClone', 'Rare_deleterious_LeukemiaGene_binom_LargeClone', 'Rare_deleterious_TOPMed_CHIPvar_binom_LargeClone','Rare_deleterious_SOMATIC_binom_LargeClone')),]


p <- ggplot(subset, aes(x = age, y = VarCount, group = Filtration, color=as.factor(Filtration)))
p <- p +geom_smooth()#+ylim(0,0.4)#method = lm, formula = y ~ splines::bs(x, 3))+ylim(0,0.4)
p <- p +facet_wrap(~Filtration, scales="free",ncol=4)+  xlab("Age") + ylab("Prevalence")
p <- p + theme(
	strip.background = element_blank(),
	strip.text.x = element_text(size = 14),
	strip.text.y = element_text(size = 14),
	axis.text.x = element_text(size=14, hjust=1, color = "black"),
	axis.text.y = element_text(size=14,hjust=1, color = "black"),
	axis.ticks =  element_line(colour = "black"), 
	axis.title.y= element_text(size=20),
	axis.title.x= element_text(size=20),
	legend.position="right", 
	legend.title = element_blank(),
	panel.background = element_blank(), 
	panel.border = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	panel.spacing = unit(1.0, "lines"), 
	plot.background = element_blank(), 
	plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
	axis.line.x = element_line(colour = "black"),
	axis.line.y = element_line(colour = "black")
)+theme(legend.position="none")

pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/LeukemiaVar.ageAssoc_spline.subset.pdf",sep=""), width = 18, height= 6)
print(p)
dev.off()

p <- ggplot(subset, aes(x = age, y = VarCount, group = Filtration, color=as.factor(Filtration)))
p <- p +geom_smooth()#+ylim(0,0.4)#method = lm, formula = y ~ splines::bs(x, 3))+ylim(0,0.4)
p <- p  +facet_wrap(~Filtration, scales="free") + xlab("Age") + ylab("Variant Count")
p <- p + theme(
	strip.background = element_blank(),
	strip.text.x = element_text(size = 14),
	strip.text.y = element_text(size = 14),
	axis.text.x = element_text(size=10, hjust=1, color = "black"),
	axis.text.y = element_text(size=12,hjust=1, color = "black"),
	axis.ticks =  element_line(colour = "black"), 
	axis.title.y= element_text(size=14),
	axis.title.x= element_text(size=14),
	legend.position="right", 
	legend.title = element_blank(),
	panel.background = element_blank(), 
	panel.border = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	panel.spacing = unit(1.0, "lines"), 
	plot.background = element_blank(), 
	plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
	axis.line.x = element_line(colour = "black"),
	axis.line.y = element_line(colour = "black")
)

pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/LeukemiaVar.ageAssoc_splineLinReg.pdf",sep=""), width = 10, height= 5)
print(p)
dev.off()

p <- ggplot(pheno9, aes(x = factor(hasCHIP), y = SOMATIC))
p <- p +geom_violin(width=1.2, fill="blue")#+ylim(0,0.4)#method = lm, formula = y ~ splines::bs(x, 3))+ylim(0,0.4)
p <- p  + xlab("CHIP") + ylab("Rare Deleterious\nKnown Somatic Variants (N)")
p <- p + theme(
	strip.background = element_blank(),
	strip.text.x = element_text(size = 14),
	strip.text.y = element_text(size = 14),
	axis.text.x = element_text(size=10, hjust=1, color = "black"),
	axis.text.y = element_text(size=12,hjust=1, color = "black"),
	axis.ticks =  element_line(colour = "black"), 
	axis.title.y= element_text(size=14),
	axis.title.x= element_text(size=14),
	legend.position="none", 
	legend.title = element_blank(),
	panel.background = element_blank(), 
	panel.border = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	panel.spacing = unit(1.0, "lines"), 
	plot.background = element_blank(), 
	plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
	axis.line.x = element_line(colour = "black"),
	axis.line.y = element_line(colour = "black")
)

pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/SOMATICVar.CHIPAssoc_boxplot.pdf",sep=""), width = 5, height= 5)
print(p)
dev.off()

subset = counts_and_phenos_LONG[which(counts_and_phenos_LONG$Filtration %in%c('Rare_deleterious', 'Rare_deleterious_LeukemiaGene', 'TOPMed_CHIPvar', 'Rare_deleterious_SOMATIC')),]
p <- ggplot(subset, aes(x = AF, y = VarCount))
#p <- p +geom_violin(width=1.2, fill="blue")#+ylim(0,0.4)#method = lm, formula = y ~ splines::bs(x, 3))+ylim(0,0.4)
p <- p  +geom_smooth()+ xlab("CHIP VAF") + ylab("Prevalence")+facet_wrap(~Filtration, scales="free",ncol=4)
p <- p + theme(
	strip.background = element_blank(),
	strip.text.x = element_text(size = 18),
	strip.text.y = element_text(size = 18),
	axis.text.x = element_text(size=18, hjust=1, color = "black"),
	axis.text.y = element_text(size=18,hjust=1, color = "black"),
	axis.ticks =  element_line(colour = "black"), 
	axis.title.y= element_text(size=20),
	axis.title.x= element_text(size=20),
	legend.position="none", 
	legend.title = element_blank(),
	panel.background = element_blank(), 
	panel.border = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	panel.spacing = unit(1.0, "lines"), 
	plot.background = element_blank(), 
	plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
	axis.line.x = element_line(colour = "black"),
	axis.line.y = element_line(colour = "black")
)
pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/VarCount.CHIPAssoc.geomSmooth.subset.pdf",sep=""), width = 18, height= 6)
print(p)
dev.off()

counts_and_phenos_LONG$hasCHIP = ifelse(counts_and_phenos_LONG$AF>0,1,0)
p <- ggplot(counts_and_phenos_LONG, aes(x = factor(hasCHIP), y = VarCount))
p <- p +geom_violin(width=1.2, fill="blue")#+ylim(0,0.4)#method = lm, formula = y ~ splines::bs(x, 3))+ylim(0,0.4)
p <- p  + xlab("CHIP") + ylab("Variants (N)")+facet_wrap(~Filtration, scales="free",ncol=4)
p <- p + theme(
	strip.background = element_blank(),
	strip.text.x = element_text(size = 14),
	strip.text.y = element_text(size = 14),
	axis.text.x = element_text(size=10, hjust=1, color = "black"),
	axis.text.y = element_text(size=12,hjust=1, color = "black"),
	axis.ticks =  element_line(colour = "black"), 
	axis.title.y= element_text(size=14),
	axis.title.x= element_text(size=14),
	legend.position="none", 
	legend.title = element_blank(),
	panel.background = element_blank(), 
	panel.border = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	panel.spacing = unit(1.0, "lines"), 
	plot.background = element_blank(), 
	plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
	axis.line.x = element_line(colour = "black"),
	axis.line.y = element_line(colour = "black")
)
pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/VarCount.CHIPAssoc.violin.pdf",sep=""), width = 18, height= 18)
print(p)
dev.off()

p <- ggplot(pheno8, aes(x = factor(hasCHIP), y = LeukemiaGene))
p <- p +geom_violin(width=1.2, fill="blue")#+ylim(0,0.4)#method = lm, formula = y ~ splines::bs(x, 3))+ylim(0,0.4)
p <- p  + xlab("CHIP") + ylab("Rare Deleterious\nin Known Leukemic Genes (N)")
p <- p + theme(
	strip.background = element_blank(),
	strip.text.x = element_text(size = 14),
	strip.text.y = element_text(size = 14),
	axis.text.x = element_text(size=10, hjust=1, color = "black"),
	axis.text.y = element_text(size=12,hjust=1, color = "black"),
	axis.ticks =  element_line(colour = "black"), 
	axis.title.y= element_text(size=14),
	axis.title.x= element_text(size=14),
	legend.position="none", 
	legend.title = element_blank(),
	panel.background = element_blank(), 
	panel.border = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	panel.spacing = unit(1.0, "lines"), 
	plot.background = element_blank(), 
	plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
	axis.line.x = element_line(colour = "black"),
	axis.line.y = element_line(colour = "black")
)

pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/LeukemicGeneVar.CHIPAssoc_boxplot.pdf",sep=""), width = 5, height= 5)
print(p)
dev.off()

counts_and_phenos_LONG$Filtrationv2 = ifelse(counts_and_phenos_LONG$Filtration == "hasCHIP", "CHIP", 
									ifelse(counts_and_phenos_LONG$Filtration == "LeukemiaGene", "Rare Deleterious,\n Leukemia Gene",
									ifelse(counts_and_phenos_LONG$Filtration == "SOMATIC", "Rare Deleterious,\n Known Somatic",
									ifelse(counts_and_phenos_LONG$Filtration == "annotated_overlap", "Rare Deleterious",
									ifelse(counts_and_phenos_LONG$Filtration == "num_FilterMutect", "QCed with Filter Mutect", NA))))) 	
counts_and_phenos_LONG$Filtrationv2 = factor(counts_and_phenos_LONG$Filtrationv2, levels=c("CHIP", "Rare Deleterious,\n Leukemia Gene","Rare Deleterious,\n Known Somatic", "Rare Deleterious", "QCed with Filter Mutect"))

p <- ggplot(counts_and_phenos_LONG[-which(counts_and_phenos_LONG$Filtration == "TOPMed_CHIPVar"),], aes(x = age, y = VarCount, group = Filtrationv2, color=as.factor(Filtrationv2)))
p <- p +geom_smooth()#+ylim(0,0.4)#method = lm, formula = y ~ splines::bs(x, 3))+ylim(0,0.4)
p <- p  +facet_wrap(~Filtrationv2, scales="free") + xlab("Age") + ylab("Variant Count")
p <- p + theme(
	strip.background = element_blank(),
	strip.text.x = element_text(size = 14),
	strip.text.y = element_text(size = 14),
	axis.text.x = element_text(size=10, hjust=1, color = "black"),
	axis.text.y = element_text(size=12,hjust=1, color = "black"),
	axis.ticks =  element_line(colour = "black"), 
	axis.title.y= element_text(size=14),
	axis.title.x= element_text(size=14),
	legend.position="none", 
	legend.title = element_blank(),
	panel.background = element_blank(), 
	panel.border = element_blank(), 
	panel.grid.major = element_blank(), 
	panel.grid.minor = element_blank(), 
	panel.spacing = unit(1.0, "lines"), 
	plot.background = element_blank(), 
	plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
	axis.line.x = element_line(colour = "black"),
	axis.line.y = element_line(colour = "black")
)

pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/LeukemiaVar.ageAssoc_splineLinReg_ALL.pdf",sep=""), width = 10, height= 5)
print(p)
dev.off()

pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/LeukemiaVar.ageAssoc.pdf",sep=""), width = 5, height= 5)
plot = ggplot(data = counts_and_phenos_LONG[-which(counts_and_phenos_LONG$Filtration == "TOPMed_CHIPVar"),], aes(x=VarCount,color= factor(Filtrationv2)))
plot=plot + geom_histogram(aes(VarCount))+xlab("Variant Count")+ylab("Percent of Participants (%)")+ scale_y_continuous(labels = scales::percent_format())
plot=plot +  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +facet_wrap(~Filtration,scales = "free") 
plot = plot  +theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12),
    panel.border = element_blank(), 
    strip.text.y = element_text(size = 12),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12,hjust=1),
    axis.ticks =  element_line(colour = "black"), 
    axis.title.x= element_text(size=12),
    axis.title.y= element_text(size=12),
    panel.background = element_blank(), 
    legend.text= element_text(size=12),
    legend.title=element_text(size=12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.margin = unit(1.0, "lines"), 
    legend.position="bottom",
    plot.background = element_blank(), 
    plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
    axis.line = element_line(colour = "black"))
p1 = plot
p1
dev.off()

pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/LeukemiaVar.pdf",sep=""), width = 10, height= 10)
plot = ggplot(data = counts_and_phenos_LONG, aes(x=age, y=VarCount, group=Filtration))
plot=plot + geom_point(aes(x=age, y=VarCount))+geom_smooth(method=lm) +xlab("AGE")+ylab("Variants (N)")
plot=plot +  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +facet_wrap(~Filtration,scales = "free") 
plot = plot  +theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12,hjust=1),
    axis.ticks =  element_line(colour = "black"), 
    axis.title.x= element_text(size=12),
    axis.title.y= element_text(size=12),
    panel.background = element_blank(), 
    panel.border = element_blank(), 
    legend.text= element_text(size=12),
    legend.title=element_text(size=12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.margin = unit(1.0, "lines"), 
    legend.position="bottom",
    plot.background = element_blank(), 
    plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
    axis.line = element_line(colour = "black"))
p1 = plot
p1
dev.off()

pdf(paste("/medpop/esp2/mzekavat/CHIP/CHUD/data/LeukemiaVar.hist.pdf",sep=""), width = 10, height= 10)
plot = ggplot(data = counts_and_phenos_LONG[which], aes(VarCount, colour=Filtration))
plot=plot + geom_histogram(aes(VarCount))+xlab("Variants (N)")+ylab("Count")
plot=plot +  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1)) +facet_wrap(~Filtration,scales = "free") 
plot = plot  +theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12,hjust=1),
    axis.ticks =  element_line(colour = "black"), 
    axis.title.x= element_text(size=12),
    axis.title.y= element_text(size=12),
    panel.background = element_blank(), 
    panel.border = element_blank(), 
    legend.text= element_text(size=12),
    legend.title=element_text(size=12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.margin = unit(1.0, "lines"), 
    legend.position="bottom",
    plot.background = element_blank(), 
    plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
    axis.line = element_line(colour = "black"))
p1 = plot
p1
dev.off()

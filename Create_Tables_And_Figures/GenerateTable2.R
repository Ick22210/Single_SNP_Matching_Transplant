#Code to Generate Table 2

mafs = c(0.05,0.1,0.2,0.3,0.4,0.5)

#for IBS mismatch score
#S = 0 if genotypes are {(aa,aa),(Aa,Aa),(AA,AA)}
#S = 1 if genotypes are {(Aa,aa),(aa,Aa),(AA,Aa),(Aa,AA)}
#S = 2 if genotypes are {(AA,aa),(aa,AA)}

#aa, aa
aa = mafs*mafs
Aa = 2*mafs*(1-mafs)
AA = (1-mafs)*(1-mafs)

pSIBS0 = (aa*aa + Aa*Aa + AA*AA)
pSIBS1 = (Aa*aa + aa*Aa + AA*Aa + Aa*AA)
pSIBS2 = (AA*aa + aa*AA)

#check that they add to 1
#pSIBS0 + pSIBS1 + pSIBS2

#For Incomp Score
#S = 0 if genotypes are {(aa,aa),(Aa,Aa),(AA,AA)}
#S = 1 if genotypes are {(Aa,aa),(aa,Aa),(AA,Aa),(Aa,AA),(AA,aa),(aa,AA)}

pSIncomp0 = (aa*aa + Aa*Aa + AA*AA)
pSIncomp1 = (Aa*aa + aa*Aa + AA*Aa + Aa*AA) + (AA*aa + aa*AA)

#check that they add to 1
#pSIncomp0 + pSIncomp1

#For MM Score
#S = 0 if genotypes are {(aa,aa),(aa,Aa),(Aa,Aa),(AA,Aa),(AA,AA)}
#S=1 if genotypes are {(Aa,aa), (Aa,AA)}
#S = 2 if genotypes are {(AA,aa),(aa,AA)}

pSMM0 = (aa*aa + aa*Aa + Aa*Aa + AA*Aa + AA*AA)
pSMM1 = (Aa*aa + Aa*AA)
pSMM2 = (AA*aa + aa*AA)

#check that they add to 1
#pSMM0 + pSMM1 + pSMM2

#merge into table
table2 = rbind(pSIBS0, pSIBS1, pSIBS2, pSIncomp0, pSIncomp1, pSMM0, pSMM1, pSMM2)
colnames(table2) = mafs

round(table2, digits = 3)

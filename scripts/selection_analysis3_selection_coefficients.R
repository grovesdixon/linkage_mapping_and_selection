#CALCULATE SELECTION COEFFICIENTS

#LOAD THE DATA OUTPUT FROM selection_analysis_setup.R
path.2.working.dir = "~/git_Repositories/linkage_mapping_and_selection/data_files/"
setwd(path.2.working.dir)
load('LinkageMappingSetup.R')

#Previous quantifications of selection were based on counts of alleles
#Here we quantify selection in terms of changes in genotype frequencies

#UPLOAD THE PARENT GENOTYPES (These are the concensus genotypes for the triplicately genotyped parent samples)
parents = read.table('concensus_parent_genotypes.txt', header = T)
colnames(parents) = c('variant', 'A', 'C')
head(parents)

#UPLOAD THE REFERENCE ALLELES (These were output as part of the genotyping pipeline)
reference.alleles = read.table('ReferenceAlleles.txt', header = T)

#REASSIGN THE VARIANT NAMES AND MAP LOCATIONS FROM THE PRIMARY DATAFRAME
variants = dat$variant
cM = dat$cM

#FUNCTION get.frequencies()
#This function pulls allele frequency data from the original counts dataframe
get.frequencies = function(dat, cTotCol, hTotCol, c.ref, h.ref, c.alt, h.alt){
  PC = c()
  PH = c()
  QC = c()
  QH = c()
  Sme = c()
  Smish = c()
  negatively.selected = c()
  underSelection = c()
  for (i in seq(1:length(dat[,cTotCol]))){
    c.tot = dat[i, cTotCol]
    h.tot = dat[i, hTotCol]
    c.a1 = dat[i, c.ref] #get reference allele count
    h.a1 = dat[i, h.ref]
    c.a2 = dat[i, c.alt]
    h.a2 = dat[i, h.alt]
    heat = c(h.a1, h.a2)
    control = c(c.a1, c.a2)
    pc = c.a1/c.tot
    ph = h.a1/h.tot
    qc = c.a2/c.tot
    qh = h.a2/h.tot
    PC = append(PC, pc)
    PH = append(PH, ph)
    QC = append(QC, qc)
    QH = append(QH, qh)
  }
  freq.dat = data.frame(PC, PH, QC, QH, cM, variants)
  colnames(freq.dat) = c('pc', 'ph', 'qc', 'qh', 'cM','variant')
  freq.dat2 = merge(parents,merge(freq.dat, reference.alleles, by = "variant"), by = 'variant')
  return(freq.dat2[with(freq.dat2, order(cM)),])
}

#FUNCTION get.s
#Calculates a selection coefficient for each variant based on the 
#change in allele frequency between heat and control replicates.
#The input is a dataframe output from get.frequencies()

#Brief Explanation:
# Since in nearly all cases the parental genotypes were of the type AA:Aa, 
#the selection process in the F1 can be viewed as a result of competition between AA homozygote and Aa heterozygote.
#In such case we calculated selection coefficients using the following equation:
# s = (p - p') / (p(1 - p'))
get.s = function(df){
  S = c()
  Sh = c()
  #We don't calculate selection coefficients for comparisons where both parents were heterozygous
  #So enter NA for those variants and skip to next
  for (i in seq(1, nrow(df))){
    A = as.character(df[i,'A'])
    C = as.character(df[i, 'C'])
    if(A == C){
      sh = NA
      Sh = append(Sh, sh)
      next
    }
    p.allele = 'ref'
    ##assign the variables
    p.c = df[i,'pc'] ; p.h = df[i,'ph'] ; qc = df[i, 'qc'] ; qh = df[i, 'qh']
    ##assign p as the shared allele between parents (this will always be the bigger allele frequency)
    pc = max(c(p.c, qc))
    ph = max(c(p.h, qh))
    qc = 1 - pc ##frequency of unshared allele pre selection
    qh = 1 - ph #frequency of unshared allele post selection
    delt.p = ph - pc ##get the change in allele frequency of the shared allele
    if (delt.p < 0) {##selection was against shared allele (negative selection on AA)
      p0 = 1 - 2*qc ##set frequency of AA genotype in terms of the unshared allele frequency
      p1 = 1 - 2*qh ##frequency of AA post selection
      sh = (p0 - p1) / (p0 *(1 - p1))
    } else{
      p0 = 2*qc ##set frequency of AT genotype in terms of the unshared allele frequency
      p1 = 2*qh ##frequency of AT post selection
      sh = (p0 - p1) / (p0 *(1 - p1))
    }##selection was for shared allele (negative selection on AT)
    Sh = append(Sh, sh)
  }
  print(max(na.omit(Sh)))
  return(Sh)
}

#FUNCTION plot.sel()
#Plots the selection coefficients across the linkage groups to make sure they align with other signals
plot.sel = function(df1, SPAN, YLIM, transform){
  df <- na.omit(df1)
  loess_fit <- loess(s ~ cM, df, span = SPAN, se = T)
  loess.dat = data.frame(df$cM, predict(loess_fit))
  colnames(loess.dat) = c('cM', 'loess')
  plot(s~cM, data = df, ylim = YLIM)
  lines(loess.dat$cM, loess.dat$loess, col = 'red', lwd = 4)
  loess.dat$loess = loess.dat$loess * transform
  return(loess.dat)
}

#RUN THE FUNCTIONS TO LOOK AT SELECTION COEFFICIENTS
ac1.s = get.frequencies(dat, "AC1c_tot", "AC1h_tot", "AC1c_refCount", "AC1h_refCount", "AC1c_altCount", "AC1h_altCount")
ac1.s$s = get.s(ac1.s)
ac1.s = na.omit(ac1.s)
ac2.s = get.frequencies(dat, "AC2c_tot", "AC2h_tot", "AC2c_refCount", "AC2h_refCount", "AC2c_altCount", "AC2h_altCount")
ac2.s$s = get.s(ac2.s)
ac2.s = na.omit(ac2.s)
ca1.s = get.frequencies(dat, "CA1c_tot", "CA1h_tot", "CA1c_refCount", "CA1h_refCount", "CA1c_altCount", "CA1h_altCount")
ca1.s$s = get.s(ca1.s)
ca1.s = na.omit(ca1.s)
ca2.s = get.frequencies(dat, "CA2c_tot", "CA2h_tot", "CA2c_refCount", "CA2h_refCount", "CA2c_altCount", "CA2h_altCount")
ca2.s$s = get.s(ca2.s)
ca2.s = na.omit(ca2.s)
YLIM = c(0, 1)
span = 0.03
transform = 1
quartz()
par(mfrow = c(4,1))
ac1.line = plot.sel(ac1.s, span, YLIM, transform)
ac2.line = plot.sel(ac2.s, span, YLIM, transform)
ca1.line = plot.sel(ca1.s, span, YLIM, transform)
ca2.line= plot.sel(ca2.s, span, YLIM, transform)

#OUTPUT THE WORKSPACE
save.image("selection_analysis3_output.R")


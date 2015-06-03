#LOOK AT CHANGE IN ALLELE FREQUENCY FOR PARTICULAR CHROMOSOMES ACROSS GENOME

#LOAD THE DATA OUTPUT FROM selection_analysis3_selection_coefficients.R
path.2.working.dir = "~/git_Repositories/linkage_mapping_and_selection/data_files/"
setwd(path.2.working.dir)
load('selection_analysis3_output.R')

#HERE WE PULL PHASING DATA OUTPUT FROM JOIN MAP AND OUTPUT IT ALONG WITH ALLELE FREQUENCY CHANGES SO WE CAN SEE 
#WHICH PARENTAL CHROMOSOMES WERE UNDER SELECTION AT OUR INDICATED REGIONS

#SWITCH TO THE PHASING DATA DIRECTORY
setwd(paste(path.2.working.dir, "/phasing_data/", sep = ""))

#MERGE THE PARENT GENOTYPES AND THE REFERENCE ALLELES INTO A SINGLE DATAFRAME
head(parents)
head(reference.alleles)
p.dat = merge(parents, reference.alleles, by = 'variant')

#HERE WE RETURN TO THE MEAN CHANGES IN ALLELE FREQUENCIES ACROSS REPLICATES
#ADD THE VARIANT NAMES AND THE MAP LOCATIONS TO THESE DATAFRAMES OF THE CHANGES IN ALLELE FREQUENCY
ac.delt.ps$variant = VARIANT
ca.delt.ps$variant = VARIANT
ac.delt.ps$cM = CM
ca.delt.ps$cM = CM


#ALSO UPDATE THE CONTROL COMPARISON
#here the delta p values were calculated as CAcontrol1 - CAcontrol2 and ACcontrol1 - ACcontrol2
con.delt.ps$variant = VARIANT
con.delt.ps$cM = CM

#BUILD A CUMULATIVE PHASING DATAFRAME FROM THE PHASING DATA FOR EACH LINKAGE GROUP
for (i in lgs){
  phase.file = paste('g', i, '.loc.phasing.txt', sep = "")
  pdat = read.table(phase.file, header = T)
  pdat$lg = i
  endpoint = endpoints[i]
  pdat$cum.pos = pdat$position + endpoint
  pdat2 = merge(pdat, p.dat, by = 'variant', all.x = TRUE)
  if (i == 1){
    cum.dat = pdat2
  }
  if (i > 1){
    cum.dat = rbind(cum.dat, pdat2)
  }
}
nrow(cum.dat)
head(cum.dat)
tail(cum.dat)

#####################################################################################
#OUTPUT THE PHASING DATA FOR PUBLICATION
#this ouput is meant to be reformatted with reformat_phasing.py
#you can then pull genomic sequencing surrounding each rad tag with pull_seq_by_RAD_tag.py
z = cum.dat[with(cum.dat, order(cum.pos)),]
head(z)
write.table(z, "compiled_phasing_data.txt", quote = F, row.names = F, sep = "\t")
#####################################################################################

#FUNCTION ref.single
#Takes vectors of the 'A' and 'C' parent genotypes and the reference allele and returns TRUE or FALSE
#for whether the reference allele is only found in one copy (example: 'A' = A/T; 'C' = T/T; ref = 'A; would return TRUE)
ref.single = function(Ageno, Cgeno, ref){
  results = c()
  for (i in 1:length(ref)){
    a = Ageno[i]
    c = Cgeno[i]
    r = ref[i]
    v = c(substr(a, 1, 1), substr(a, 3, 3), substr(c, 1, 1), substr(c, 3, 3))
    sub = v[v == r]
    if (length(sub) == 1){
      results = append(results, 'TRUE')
    }
    if (length(sub) > 1){
      results = append(results, 'FALSE')
    }
  }
  return(results)
}
#FUNCTION chrom.dp
#Returns the change in allele frequency for particular chromosomes
chrom.dp = function(dat, chromosome, ref, diplotype, ref.single, delta.p){
  #ARGUMENTS:
  #dat = the dataframe
  #chromosomes = the column with the chromosomes phasing data
  #ref = the reference allele column
  #diplotype = the column with the parent's diploid genotypes
  #ref.single = column that says whether the reference allele is single copy among parents
  #delta.p = the column with the change in allele frequency for the reference allele
  cdp = c()
  for (i in 1:nrow(dat)){
    map.allele = dat[i, chromosome]
    ref.allele = dat[i, ref]
    geno = dat[i, diplotype]
    first.allele = substr(geno, 1, 1)
    rs = dat[i, ref.single]
    dp = dat[i, delta.p]
    if (map.allele == "n"){
      if (rs == 'TRUE'){
        cdp = append(cdp, -1*dp)
        next
      }
      if (rs == 'FALSE'){
        cdp = append(cdp, dp)
        next
      }
    }
    if (map.allele == "p"){
      if (rs == 'FALSE'){
        cdp = append(cdp, -1*dp)
        next
      }
      if (rs == 'TRUE'){
        cdp = append(cdp, dp)
        next
      }
    }
    if (map.allele == "m"){
      if (rs == 'FALSE'){
        cdp = append(cdp, -1*dp)
        next
      }
      if (rs == 'TRUE'){
        cdp = append(cdp, dp)
        next
      }
    }
    if (map.allele == "l"){
      if (rs == 'TRUE'){
        cdp = append(cdp, -1*dp)
        next
      }
      if (rs == 'FALSE'){
        cdp = append(cdp, dp)
        next
      }
    }
    if (map.allele == "h"){
      if (ref.allele == first.allele){
        cdp = append(cdp, dp)
      }
      if (ref.allele != first.allele){
        cdp = append(cdp, -1*dp)
      }
    }
    if (map.allele == "k"){
      if (ref.allele == first.allele){
        cdp = append(cdp, -1*dp)
      }
      if (ref.allele != first.allele){
        cdp = append(cdp, dp)
      }
    }
  }
  return(cdp)
}

#SET UP MERGED DATASETS WITH THE FORMATTED PHASING DATA AND THE CHANGES IN ALLELE FREQUENCY
ac.dat = merge(cum.dat, ac.delt.ps, by = 'variant')
ca.dat = merge(cum.dat, ca.delt.ps, by = 'variant')
control.dat = merge(cum.dat, con.delt.ps, by = 'variant')

#USE ref.single() TO SET UP A COLUMN SAYING WHETHER THE REFERENCE ALLELE IS SINGLE COPY AMONG THE PARENTS
ac.dat$ref.single = ref.single(ac.dat$A, ac.dat$C, ac.dat$ref)
ca.dat$ref.single = ref.single(ca.dat$A, ca.dat$C, ca.dat$ref)
control.dat$ref.single = ref.single(control.dat$A, control.dat$C, control.dat$ref)

#COPY THE CONTROL DATAFRAME INTO TWO SEPARATE ONES, ONE FOR EACH CROSS
#this is because it doesn't make sense to use the average change in allele frequency
AC.control.dat = control.dat
CA.control.dat = control.dat

#CALCULATE THE ALLELE FREQUENCY CHANGES FOR HETEROZYGOUS LOCI ON EACH CHROMOSOME
#for ac replicates
ac.dat$AC1dp = chrom.dp(ac.dat, 'A.c1', 'ref', 'A', 'ref.single', 'mean')
ac.dat$AC2dp = chrom.dp(ac.dat, 'A.c2', 'ref', 'A', 'ref.single', 'mean')
ac.dat$CC1dp = chrom.dp(ac.dat, 'C.c1', 'ref', 'C', 'ref.single', 'mean')
ac.dat$CC2dp = chrom.dp(ac.dat, 'C.c2', 'ref', 'C', 'ref.single', 'mean')
#for ca replicates
ca.dat$AC1dp = chrom.dp(ca.dat, 'A.c1', 'ref', 'A', 'ref.single', 'mean')
ca.dat$AC2dp = chrom.dp(ca.dat, 'A.c2', 'ref', 'A', 'ref.single', 'mean')
ca.dat$CC1dp = chrom.dp(ca.dat, 'C.c1', 'ref', 'C', 'ref.single', 'mean')
ca.dat$CC2dp = chrom.dp(ca.dat, 'C.c2', 'ref', 'C', 'ref.single', 'mean')
#for the ac controls
AC.control.dat$AC1dp = chrom.dp(AC.control.dat, 'A.c1', 'ref', 'A', 'ref.single', 'AC.Controls.delta.p')
AC.control.dat$AC2dp = chrom.dp(AC.control.dat, 'A.c2', 'ref', 'A', 'ref.single', 'AC.Controls.delta.p')
#for the ca controls
CA.control.dat$CA1dp = chrom.dp(CA.control.dat, 'C.c1', 'ref', 'C', 'ref.single', 'CA.Controls.delta.p')
CA.control.dat$CA2dp = chrom.dp(CA.control.dat, 'C.c2', 'ref', 'C', 'ref.single', 'CA.Controls.delta.p')

#MAKE TWO SUBSET DATAFRAMES OF ONLY THE HETEROZYGOUS VARIANTS FOR EACH PARENT
# #for the ac replicates (parent A is mom)
ac.a.poly = ac.dat[ac.dat$crossType != 'nnxnp',] #include variants except when parent A is homozygous
ac.c.poly = ac.dat[ac.dat$crossType != 'lmxll',] #include variants except when parent C is homozygous
ac.a.poly = ac.a.poly[with(ac.a.poly, order(cM)),]
ac.c.poly = ac.c.poly[with(ac.c.poly, order(cM)),]
#for the ca replicates (parent C is mom)
#(note that for linkage mapping all markers were noted as if parent A were mother, so nnxnp indicates that C was heterozygous)
ca.c.poly = ca.dat[ca.dat$crossType != 'lmxll',] #include variants except when parent C is homozygous
ca.a.poly = ca.dat[ca.dat$crossType != 'nnxnp',] #include variants except when parent A is homozygous
ca.c.poly = ca.c.poly[with(ca.c.poly, order(cM)),]
ca.a.poly = ca.a.poly[with(ca.a.poly, order(cM)),]
#for the ac control comparisons
ac.control.poly = AC.control.dat[AC.control.dat$crossType != 'nnxnp',]
ca.control.poly = CA.control.dat[CA.control.dat$crossType != 'lmxll',]
ac.control.poly = ac.control.poly[with(ac.control.poly, order(cM)),]
ca.control.poly = ca.control.poly[with(ca.control.poly, order(cM)),]

#if not using only hets uncomment these
# ac.a.poly = ac.dat[with(ac.dat, order(cM)),]
# ac.c.poly = ac.dat[with(ac.dat, order(cM)),]
# ca.a.poly = ca.dat[with(ac.dat, order(cM)),]
# ca.c.poly = ca.dat[with(ac.dat, order(cM)),]

#FUNCTION plot.chroms
#Plots the chromosome specific allele frequency changes for a given parent
plot.chroms = function(dat, lgs, chrom1.col, chrom2.col, YLIM, MAIN, XLAB, YLAB, COLORS, CUT, ADD){
  if (ADD == FALSE){
    plot(dat[,chrom1.col]~dat$cM, type = "n", ylim = YLIM, main = MAIN, axes = F, xlab = XLAB, ylab = YLAB)
  }
  for (i in lgs){
    lg.sub = dat[dat$lg == i,]
    lines(lg.sub[,chrom1.col]~lg.sub$cM, col = COLORS[1], lwd = 1)
    lines(lg.sub[,chrom2.col]~lg.sub$cM, col = COLORS[2], lwd = 1)
  }
  abline(v = endpoints, col = 'grey', lwd = 0.5)
  sig = dat[p.adjust(dat$fisher.p, method = "BH") < CUT,]
  points(abs(sig[,chrom1.col])~sig$cM, pch = 19, col = 'red', cex = 0.5)
  #   legend(1410, 0, c("Chromosome 1", "Chromosome 2"), fill = c('red', 'blue'))
  axis(1, at = endpoints, labels = F)
  axis(2, las = 1)
  mtext(lg.names, side = 1, at = mids, line = .75)
}

#UPLOAD BOOTSTRAP BAR COORDINATES FROM selection_analysis2_chiSquare.R
setwd(path.2.working.dir)
ac.bar = read.table("ac_bars.txt")
ca.bar = read.table("ca_bars.txt")

#FUNCTION TO PLOT THE BARS 
add.bars = function(bar.dat, bar.colors){
  easy = bar.dat[bar.dat$cut == 'easy',]
  for (i in 1:nrow(easy)){
    x0 = easy$lefts[i]
    x1 = easy$rights[i]
    y = .0
    segments(x0, y, x1, y, lwd = 3, col = bar.colors[1])
  }
  hard = bar.dat[bar.dat$cut == 'hard',]
  print(hard)
  for (i in 1:nrow(hard)){
    x0 = hard$lefts[i]
    x1 = hard$rights[i]
    y = 0
    segments(x0, y, x1, y, lwd = 3, col = bar.colors[2])
  }
}

######################################################################################
############## PLOT THE MEAN ALLELE FREQUENCY AVERAGED ACROSS REPLICATES #############
######################################################################################

#PLOT FOR THE AC REPLICATES
YLIM = c(-.3, .3)
COLORS = c('cyan3', 'coral')
bar.colors = c('darkgreen', 'yellowgreen')
CUT = 0.05
quartz()
par(mfrow = c(2,1))
titleM = "AC (heat vs control) Maternal Chromosomes"
titleP = "AC (heat vs control) Paternal Chromosomes"
YLAB = expression(paste('P'['heat'], " - P"["control"], sep = ""))
plot.chroms(ac.a.poly, lgs, 'AC1dp', 'AC2dp', YLIM, titleM, "Linkage Group", YLAB, COLORS, CUT, FALSE)
plot.chroms(ac.c.poly, lgs, 'CC1dp', 'CC2dp', YLIM, titleP, "Linkage Group", YLAB, COLORS, CUT, FALSE)
add.bars(ac.bar, bar.colors)

#PLOT FOR THE CA REPLICATES
quartz()
par(mfrow = c(2,1))
titleM = "CA (heat vs control) Maternal Chromosomes"
titleP = "CA (heat vs control) Paternal Chromosomes"
plot.chroms(ca.c.poly, lgs, 'CC1dp', 'CC2dp', YLIM, titleM, "Linkage Group", YLAB, COLORS, CUT, FALSE)
plot.chroms(ca.a.poly, lgs, 'AC1dp', 'AC2dp', YLIM, titleP, "Linkage Group", YLAB, COLORS, CUT, FALSE)
add.bars(ca.bar, bar.colors)

#PLOT FOR THE CONTROLS
quartz()
par(mfrow = c(2,1))
titleM = "AC (control vs control) A Chromosomes"
titleP = "CA (control vs control) C Chromosomes"
plot.chroms(ac.control.poly, lgs, 'AC1dp', 'AC2dp', YLIM, titleM, "Linkage Group", YLAB, COLORS, CUT, FALSE)
plot.chroms(ca.control.poly, lgs, 'CA1dp', 'CA2dp', YLIM, titleP, "Linkage Group", YLAB, COLORS, CUT, FALSE)
##################################################################################################
##################################################################################################

#ABOVE WE LOOKED AT THE AVERAGES ACROSS REPLICATES
#NOW LOOK AT EACH REPLICATE INDIVIDUALLY

#CALCULATE TEH CHROMOSOME-SPECIFIC ALLELE FREQUENCY CHANGES 
#for ac replicates 1 and 2
ac.dat$AC1dp.rep1 = chrom.dp(ac.dat, 'A.c1', 'ref', 'A', 'ref.single', 'ac1.ps.delt.p')
ac.dat$AC2dp.rep1 = chrom.dp(ac.dat, 'A.c2', 'ref', 'A', 'ref.single', 'ac1.ps.delt.p')
ac.dat$AC1dp.rep2 = chrom.dp(ac.dat, 'A.c1', 'ref', 'A', 'ref.single', 'ac2.ps.delt.p')
ac.dat$AC2dp.rep2 = chrom.dp(ac.dat, 'A.c2', 'ref', 'A', 'ref.single', 'ac2.ps.delt.p')

ac.dat$CC1dp.rep1 = chrom.dp(ac.dat, 'C.c1', 'ref', 'C', 'ref.single', 'ac1.ps.delt.p')
ac.dat$CC2dp.rep1 = chrom.dp(ac.dat, 'C.c2', 'ref', 'C', 'ref.single', 'ac1.ps.delt.p')
ac.dat$CC1dp.rep2 = chrom.dp(ac.dat, 'C.c1', 'ref', 'C', 'ref.single', 'ac2.ps.delt.p')
ac.dat$CC2dp.rep2 = chrom.dp(ac.dat, 'C.c2', 'ref', 'C', 'ref.single', 'ac2.ps.delt.p')
#for ca replicates
ca.dat$AC1dp.rep1 = chrom.dp(ca.dat, 'A.c1', 'ref', 'A', 'ref.single', 'ca1.ps.delt.p')
ca.dat$AC2dp.rep1 = chrom.dp(ca.dat, 'A.c2', 'ref', 'A', 'ref.single', 'ca1.ps.delt.p')
ca.dat$AC1dp.rep2 = chrom.dp(ca.dat, 'A.c1', 'ref', 'A', 'ref.single', 'ca2.ps.delt.p')
ca.dat$AC2dp.rep2 = chrom.dp(ca.dat, 'A.c2', 'ref', 'A', 'ref.single', 'ca2.ps.delt.p')

ca.dat$CC1dp.rep1 = chrom.dp(ca.dat, 'C.c1', 'ref', 'C', 'ref.single', 'ca1.ps.delt.p')
ca.dat$CC2dp.rep1 = chrom.dp(ca.dat, 'C.c2', 'ref', 'C', 'ref.single', 'ca1.ps.delt.p')
ca.dat$CC1dp.rep2 = chrom.dp(ca.dat, 'C.c1', 'ref', 'C', 'ref.single', 'ca2.ps.delt.p')
ca.dat$CC2dp.rep2 = chrom.dp(ca.dat, 'C.c2', 'ref', 'C', 'ref.single', 'ca2.ps.delt.p')

#FILTER THE DATASETS DOWN TO ONLY THE HETEROZYGOUS LOCI AND ORDER THEM
#for the ac replicates (parent A is mom)
ac.a.poly = ac.dat[ac.dat$crossType != 'nnxnp',] #include variants except when parent A is homozygous
ac.c.poly = ac.dat[ac.dat$crossType != 'lmxll',] #include variants except when parent C is homozygous
ac.a.poly = ac.a.poly[with(ac.a.poly, order(cM)),]
ac.c.poly = ac.c.poly[with(ac.c.poly, order(cM)),]
#for the ca replicates (parent C is mom)
#(note that for linkage mapping all markers were noted as if parent A were mother, so nnxnp indicates that C was heterozygous)
ca.c.poly = ca.dat[ca.dat$crossType != 'lmxll',] #include variants except when parent C is homozygous
ca.a.poly = ca.dat[ca.dat$crossType != 'nnxnp',] #include variants except when parent A is homozygous
ca.c.poly = ca.c.poly[with(ca.c.poly, order(cM)),]
ca.a.poly = ca.a.poly[with(ca.a.poly, order(cM)),]
#for the ac control comparisons
ac.control.poly = AC.control.dat[AC.control.dat$crossType != 'nnxnp',]
ca.control.poly = CA.control.dat[CA.control.dat$crossType != 'lmxll',]
ac.control.poly = ac.control.poly[with(ac.control.poly, order(cM)),]
ca.control.poly = ca.control.poly[with(ca.control.poly, order(cM)),]

#IF NOT USING JUST HETS
# ac.a.poly = ac.dat[with(ac.dat, order(cM)),]
# ac.c.poly = ac.dat[with(ac.dat, order(cM)),]
# ca.a.poly = ca.dat[with(ac.dat, order(cM)),]
# ca.c.poly = ca.dat[with(ac.dat, order(cM)),]

#################################################################
################# PLOT WITH REPLICATES OVERLAID #################
#################################################################
#cross AC parent A
quartz()
COLORS = c('cyan3', 'coral')
bar.colors = c('darkgreen', 'yellowgreen')
par(mfrow = c(1,1))
YLIM = c(-.4, .4)
plot.chroms(ac.a.poly, lgs, 'AC1dp.rep1', 'AC2dp.rep1', YLIM, 'AC cross, A chromosomes (dam)', "Linkage Group", YLAB, COLORS, CUT, FALSE)
plot.chroms(ac.a.poly, lgs, 'AC1dp.rep2', 'AC2dp.rep2', YLIM, 'AC cross, A chromosomes (dam)', "Linkage Group", YLAB, COLORS, CUT, TRUE)
#cross AC parent C
plot.chroms(ac.c.poly, lgs, 'CC1dp.rep1', 'CC2dp.rep1', YLIM, 'AC cross, C chromosomes (sire)', "Linkage Group", YLAB, COLORS, CUT, FALSE)
plot.chroms(ac.c.poly, lgs, 'CC1dp.rep2', 'CC2dp.rep2', YLIM, 'AC cross, C chromosomes (sire)', "Linkage Group", YLAB, COLORS, CUT, TRUE)
add.bars(ac.bar, bar.colors)
#cross CA parent C
plot.chroms(ca.c.poly, lgs, 'CC1dp.rep1', 'CC2dp.rep1', YLIM, 'CA cross, C chromosomes (dam)', "Linkage Group", YLAB, COLORS, CUT, FALSE)
plot.chroms(ca.c.poly, lgs, 'CC1dp.rep2', 'CC2dp.rep2', YLIM, 'CA cross, C chromosomes (dam)', "Linkage Group", YLAB, COLORS, CUT, TRUE)
#cross CA parent A
plot.chroms(ca.a.poly, lgs, 'AC1dp.rep1', 'AC2dp.rep1', YLIM, 'CA cross, A chromosomes (sire)', "Linkage Group", YLAB, COLORS, CUT, FALSE)
plot.chroms(ca.a.poly, lgs, 'AC1dp.rep2', 'AC2dp.rep2', YLIM, 'CA cross, A chromosomes (sire)', "Linkage Group", YLAB, COLORS, CUT, TRUE)
add.bars(ca.bar, bar.colors)



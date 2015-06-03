########################################################
################# selection_analysis.R #################
########################################################
#THIS SCRIPT IS FOR STATISTICAL ANALYSIS AND FIGURE GENERATION BASED ON ALLELE 
#FREQUENCY DIFFERENCES BETWEEN HEAT SELECTED AND CONTROL LARVAL COHORTS
#Groves Dixon
#written 9-2-14
#last updated 3-3-15

####################################################
######### IMPORT AND FORMAT THE INPUT DATA #########
####################################################
path.2.working.dir = "~/git_Repositories/linkage_mapping_and_selection/data_files/"
setwd(path.2.working.dir)

#IMPORT ALLELE COUNTS
fdat = read.table('selection_analysis_Rinput.txt', header = T)
colnames(fdat)[1] = "variant"
head(fdat)
  #Allele count data are given for each replicate pair
  # 'AC' represents one cross 'CA' the other cross
  # 'AC1c' and 'AC2c' are the first and second control replicate for hte AC cross
  # 'AC1h' and 'AC2h' are the first and second heat selected replicates for hte AC cross
  # The numbers in the columns are the counts for either the reference or alternative allele calls in each replicate

#IMPORT THE MAP
map = read.table('linkage_map_Rinput.txt', header = T)
head(map)
#(the cM values are cumulative across linkage groups)
tail(map)

#MERGE THE DATAFRAMES TO PAIR VARIANTS' ALLELE COUNTS WITH MAP LOCATIONS
dat = merge(fdat, map, by = 'variant')
dat = dat[with(dat, order(cM)),]
head(dat)
length(dat$variant)

#SET UP LABELS FOR THE 14 LINKAGE GROUPS
lgs = 1:14

#GATHER THE ENDPOINTS AND MIDPOINTS OF EACH LINKAGE GROUP FOR PLOTTING PURPOSES
endpoints = c(0)
for (i in lgs){
  sub = dat[dat$lg == i,]
  end = max(sub$cM)
  endpoints = append(endpoints, end)
}
mids = c()
for (i in seq(1:(length(endpoints) - 1))){
  left = endpoints[i]
  right = endpoints[(i + 1)]
  point = (left + right) / 2
  mids = append(mids, point)
}


######################################################################
############# SET UP SOME GLOBAL VARIABLES AND FUNCTIONS #############
######################################################################
CM = dat$cM
LG = dat$lg
VARIANT = dat$variant
lg.names = 1:14
point.size = 1
setYlim = c(-.2, 8)

#FUNCTION fisher.method():
#performs fisher's method to combine a set of p values
#ARGUMENTS: ps = a vector of p values you want combined
fisher.method=function(ps) {
  ftest=-2*sum(log(ps))
  df=length(2*ps)
  pv=1-pchisq(q=ftest,df=df)
}

#FUNCTION get.combined():
#runs fisher.method() on a dataframe of p values
#ARGUMENTS: DF = a two column data frame with p values in it
get.combined = function(DF){
  combined.p.values = c()
  alt = c()
  for (i in seq(1, length(DF[,1]))){ #(apparently did not know of nrow() yet)
    x = DF[i,1]
    y = DF[i,2]
    combined.p = fisher.method(c(x,y))
    combined.p.values = append(combined.p.values, combined.p)
  }
  return(combined.p.values)
}

#FUNCTION is.odd()
# returns T/F whether x is odd
is.odd <- function(x) x %% 2 != 0


#FUNCTION plot.scan()
#Adjusts a set of p values for FDR and plots them based on map locations and adds bars below to indicate boostrap significances
#ARGUMENTS: 
  # DF = the dataframe
  # MAIN = the title you want for the plot
  # LEFT = vector of the left sides of bootstrap bars
  # Y = the Y coordinate for where to plot the bootstrap bars
  # LEN = vector of lengths for the bootstrap bars
plot.scan = function(DF, MAIN, LEFT, Y, LEN){##NOTE THIS FUNCTION IS REPEATED BELOW
  CEX = 1
  adj = p.adjust(DF$p, method = 'BH')
  #set up grey and black points alternating by LG
  colors = is.odd(DF$LG)
  colors[colors == TRUE] <- 'black'
  colors[colors == FALSE] <- 'grey'
  #make significant points red
  x = adj < 0.05
  for (i in 1:length(x)){
    if (x[i] == TRUE){
      colors[i] <- 'red'
    }
  }
  plot(-log(p, 10)~CM, data = DF, main = MAIN, pch = 19, col = colors, cex = point.size, axes = F, xlab = "Linkage Group", ylab = expression('-log'['10']*'(p)'), ylim=setYlim)
  axis(1, at = endpoints, labels = F)
  axis(2, las = 1)
  mtext(lg.names, side = 1, at = mids, line = .75)
  par(lend = 2)
  segments(LEFT, Y, (LEFT + LEN), Y, lwd = 2)
  TEXT = paste(LEN, "cM")
  labPos = LEFT + .5 * LEN
  text(labPos, y = (Y + .4), labels = TEXT)
}


#FUNCTION do.plot()
#Runs get.combined and plot.scan to build a Manhattan plot and return the modified plot data
do.plot = function(pval.df, MAIN, LEFT, Y, LEN){
  p = get.combined(pval.df)
  plot.dat = data.frame(VARIANT, LG, CM, p)
  plot.scan(plot.dat, MAIN, LEFT, Y, LEN)
  rownames(plot.dat) = plot.dat$VARIANT
  return(plot.dat)
}

#FUNCTION get.dir()
#returns the change in allele frequency for a given heat and control pair
get.dir = function(dat, cTotCol, hTotCol, c.ref, h.ref, c.alt, h.alt){
  p.greater.list = c()
  p.less.list = c()
  delt.p = c()#empty set to store the signed change in allele frequency for selected samples
  peak.ps = c()
  use.greater = c()
  for (i in seq(1:length(dat[,cTotCol]))){
    c.tot = dat[i, cTotCol]
    h.tot = dat[i, hTotCol]
    c.a1 = dat[i, c.ref] #get reference allele count
    h.a1 = dat[i, h.ref]
    c.a2 = dat[i, c.alt]
    h.a2 = dat[i, h.alt]
    heat = c(h.a1, h.a2)
    control = c(c.a1, c.a2)
    dp = h.a1/h.tot - c.a1/c.tot
    delt.p = append(delt.p, dp)
  }
  p.dat = data.frame(delt.p)
  return(p.dat)
}

################################################################################################
### GENERATE THE SET OF 'ALTERNATIVES' FOR EACH LOCUS BASED ON CHANGE IN ALLELE FREQUENCY ######
################################################################################################

#Because we are using Fisher's method to combine p-values for replicates
#we need to run 1-tailed tests for significance of allele frequency differences
#otherwise replicates that gave significant but contradictory results (ie significant in opposite tails)
# would return an inappropriately strong combined p value. (Downstream we multiply the p-values by 2 to account for 
#fact that we had no  a priori expectation of direction)
#so we need an alternative hypothesis for each 1-tailed test.
#We base this off of the mean change in frequency across replicates
#below we do this for each set of comparisons we make

##GENERATE A VECTOR WITH 'ALTERNATIVES' TO USE FOR THE 1-TAILED FISHER'S TEST FOR EACH LOCUS
##get change in allele frequency for the two 'AC' replicate pairs
ac1.ps = get.dir(dat, "AC1c_tot", "AC1h_tot", "AC1c_refCount", "AC1h_refCount", "AC1c_altCount", "AC1h_altCount")
ac2.ps = get.dir(dat, "AC2c_tot", "AC2h_tot", "AC2c_refCount", "AC2h_refCount", "AC2c_altCount", "AC2h_altCount")

#assemble into a dataframe
ac.delt.ps = data.frame(ac1.ps$delt.p, ac2.ps$delt.p) ##assemble into a dataframe

#get the mean change in allele frequency across replicates
ac.delt.ps$mean = apply(ac.delt.ps, 1, mean) 

#get T/F values for whether change in allele frequency is positive or negative 
ac.use.greater = ac.delt.ps$mean <= 0

#use the T/F values to generate a vector of alternative hypothesis values ('less' or 'greater')
ac.alternative = ac.use.greater
ac.alternative[ac.alternative == TRUE] <- 'greater'
ac.alternative[ac.alternative == FALSE] <- 'less'
###Now we have our set of alternative hypotheses for the 'AC' replicates

###REPEAT FOR THE CA REPLICATES
ca1.ps = get.dir(dat, "CA1c_tot", "CA1h_tot", "CA1c_refCount", "CA1h_refCount", "CA1c_altCount", "CA1h_altCount")
ca2.ps = get.dir(dat, "CA2c_tot", "CA2h_tot", "CA2c_refCount", "CA2h_refCount", "CA2c_altCount", "CA2h_altCount")
ca.delt.ps = data.frame(ca1.ps$delt.p, ca2.ps$delt.p)
ca.delt.ps$mean = apply(ca.delt.ps, 1, mean)
ca.use.greater = ca.delt.ps$mean <= 0
ca.alternative = ca.use.greater
ca.alternative[ca.alternative == TRUE] <- 'greater'
ca.alternative[ca.alternative == FALSE] <- 'less'

###REPEAT FOR CONTROL COMPARISONS COMPARING AC CONTROL REPLICATES WITH CA CONTROL REPLICATES
con1.ps = get.dir(dat, "CA1c_tot", "AC1c_tot", "CA1c_refCount", "AC1c_refCount", "CA1c_altCount", "AC1c_altCount")
con2.ps = get.dir(dat, "CA2c_tot", "AC2c_tot", "CA2c_refCount", "AC2c_refCount", "CA2c_altCount", "AC2c_altCount")
con.delt.ps = data.frame(con1.ps$delt.p, con2.ps$delt.p)
con.delt.ps$mean = apply(con.delt.ps, 1, mean) 
con.use.greater = con.delt.ps$mean <= 0 
con.alternative = con.use.greater 
con.alternative[con.alternative == TRUE] <- 'greater' 
con.alternative[con.alternative == FALSE] <- 'less'

##REPEAT COMPARING CONTROL REPLICATES TO CONTROL REPLICATES WITHIN CROSSES
con1.ps = get.dir(dat, "CA1c_tot", "CA2c_tot", "CA1c_refCount", "CA2c_refCount", "CA1c_altCount", "CA2c_altCount")
con2.ps = get.dir(dat, "AC1c_tot", "AC2c_tot", "AC1c_refCount", "AC2c_refCount", "AC1c_altCount", "AC2c_altCount")
con.delt.ps = data.frame(con1.ps$delt.p, con2.ps$delt.p)
head(con.delt.ps)
#reassign the column names so it's clear what was compared downstream
colnames(con.delt.ps) = c("CA.Controls.delta.p", "AC.Controls.delta.p")
con.delt.ps$mean = apply(con.delt.ps, 1, mean) 
con.use.greater = con.delt.ps$mean <= 0
con.alternative[con.alternative == TRUE] <- 'greater'
con.alternative[con.alternative == FALSE] <- 'less'

################################################################################################

#OUTPUT THE R DATASET TO RUN OTHER SCRIPTS ON
save.image(file = "LinkageMappingSetup.R")


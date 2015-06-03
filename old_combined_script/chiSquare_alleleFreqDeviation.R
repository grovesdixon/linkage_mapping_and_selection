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

################################################################################################

###NOW WE CAN USE THE ALTERNATIVES TO RUN ONE-SIDED FISHER TESTS FOR EACH REPLICATE PAIR

#FUNCITON get.one.side.ps()
#Returns a set of one sided p-values for a given replicate pair
#ARGUMENTS: dat, cTotCol, hTotCol, c.ref, h.ref, c.alt, and h.alt specify the dataframe and columns; ALT specifies the vector of alternatives to use.
get.one.side.ps = function(dat, cTotCol, hTotCol, c.ref, h.ref, c.alt, h.alt, ALT){
  p.list = c()
  for (i in seq(1:length(dat[,cTotCol]))){
    c.tot = dat[i, cTotCol]
    h.tot = dat[i, hTotCol]
    c.a1 = dat[i, c.ref] 
    h.a1 = dat[i, h.ref]
    c.a2 = dat[i, c.alt]
    h.a2 = dat[i, h.alt]
    heat = c(h.a1, h.a2)
    control = c(c.a1, c.a2)
    M = as.table(cbind(control, heat))
    rownames(M) = c("ref","alt")
    colnames(M) = c("control", "heated")
    loc = dat[i,'cM']
    var = dat[i, 'variant']
    #     print(var)
    #     print(M)
    DIR = ALT[i]
    p = fisher.test(M, alternative = DIR)$p.value
    p2 = p*2 ##multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
    p2[p2 > 1] <- 1
    p.list = append(p.list, p2)
  }
  return(p.list)
}

## NOW WE CAN PLOT THE MANHATTAN PLOTS FOR EACH REPLICATE
####### DO AC ONE-TAILED PLOT ##############
ac1.one.side = get.one.side.ps(dat, "AC1c_tot", "AC1h_tot", "AC1c_refCount", "AC1h_refCount", "AC1c_altCount", "AC1h_altCount", ac.alternative)
ac2.one.side = get.one.side.ps(dat, "AC2c_tot", "AC2h_tot", "AC2c_refCount", "AC2h_refCount", "AC2c_altCount", "AC2h_altCount", ac.alternative)
ac.ps = data.frame(ac1.one.side, ac2.one.side)
head(ac.ps)
par(mfrow = c(1,1))
ac.df = do.plot(ac.ps, 'Selection in AC Replicates', 100, 5, 100)
####### DO CA ONE-TAILED PLOT ##############
ca1.one.side = get.one.side.ps(dat, "CA1c_tot", "CA1h_tot", "CA1c_refCount", "CA1h_refCount", "CA1c_altCount", "CA1h_altCount", ca.alternative)
ca2.one.side = get.one.side.ps(dat, "CA2c_tot", "CA2h_tot", "CA2c_refCount", "CA2h_refCount", "CA2c_altCount", "CA2h_altCount", ca.alternative)
ca.ps = data.frame(ca1.one.side, ca2.one.side)
ca.df = do.plot(ca.ps, 'Selection in CA Replicates', 100, 5, 100)
####### DO CONTROL ONE-TAILED PLOT COMPARING CONTROLS BETWEEN CROSSES ##############
# con1.ps = get.one.side.ps(dat, "CA1c_tot", "AC1c_tot", "CA1c_refCount", "AC1c_refCount", "CA1c_altCount", "AC1c_altCount", con.alternative)
# con2.ps = get.one.side.ps(dat, "CA2c_tot", "AC2c_tot", "CA2c_refCount", "AC2c_refCount", "CA2c_altCount", "AC2c_altCount", con.alternative)
# con.ps = data.frame(con1.ps, con2.ps)
# con.df = do.plot(con.ps, 'Selection in Controls', 100, 5, 100)
####### DO WITHIN CROSS CONTROL ONE-TAILED PLOT ##############
con1.ps = get.one.side.ps(dat, "CA1c_tot", "CA2c_tot", "CA1c_refCount", "CA2c_refCount", "CA1c_altCount", "CA2c_altCount", con.alternative)
con2.ps = get.one.side.ps(dat, "AC1c_tot", "AC2c_tot", "AC1c_refCount", "AC2c_refCount", "AC1c_altCount", "AC2c_altCount", con.alternative)
con.ps = data.frame(con1.ps, con2.ps)
con.df = do.plot(con.ps, 'Selection in Controls', 100, 5, 100)
#----------------- DO ONE-TAILED PLOT ------------------------------------------
setYlim = c(-.2, 7.2)
par(mfrow = c(3,1))
ac.df = do.plot(ac.ps, 'Selection in AC Replicates', 100, 5, 100)
ca.df = do.plot(ca.ps, 'Selection in CA Replicates', 100, 5, 100)
con.df = do.plot(con.ps, 'Selection in Controls', 100, 5, 100)

#SAVE THE COMBINED P VALUES FOR EACH COMPARISON
ac.combined.p = ac.df$p
ca.combined.p = ca.df$p
con.combined.p = con.df$p

#THES PLOTS SHOW US THE PRIMARY PEAKS. BUT WE ALSO WANT TO LOOK FOR MORE SUBTLE SIGNALS
#WE DO THAT WITH A BOOTSTRAPPING ANALYSIS

##############################################################################
############ LOOK AT ABSOLUTE CHANGES IN ALLELE FREQUENCY ####################
##############################################################################
#These plots are intended to show that the raw changes in allele frequencies match
#with the p value peaks and are consistent accross replicates
plot.abs = function(DF, Ycolumn, MAIN){##NOTE THIS FUNCTION IS REPEATED BELOW
  CEX = 1
  #set up grey and black points alternating by LG
  colors = is.odd(LG)
  colors[colors == TRUE] <- 'black'
  colors[colors == FALSE] <- 'grey'
  plot(abs(DF[,Ycolumn])~CM, main = MAIN, pch = 19, col = colors, cex = point.size, axes = F, xlab = "Linkage Group", ylab =  'Abs(Pheat - Pcontrol)', ylim=NULL)
  axis(1, at = endpoints, labels = F)
  axis(2, las = 1)
  mtext(lg.names, side = 1, at = mids, line = .75)
  loess_fit <- loess(abs(DF[,Ycolumn])~CM, span = LOESS.SPAN, se = T)
  lines(CM, predict(loess_fit),col="red",lwd=1.5)
}
quartz()
LOESS.SPAN = .05
par(mfrow = c(4,1))
plot.abs(ac.delt.ps, 'ac1.ps.delt.p', 'AC1 changes in allele frequencies')
plot.abs(ac.delt.ps, 'ac2.ps.delt.p', 'AC2 changes in allele frequencies')
plot.abs(ca.delt.ps, 'ca1.ps.delt.p', 'CA1 changes in allele frequencies')
plot.abs(ca.delt.ps, 'ca2.ps.delt.p', 'CA2 changes in allele frequencies')
plot.abs(con.delt.ps)
##############################################################################

##################################################################################
################# DO BOOTSTRAP TO DEMONSTRATE HOW RARE PEAKS ARE #################
##################################################################################
###CREATING THE NULL DISTRIBUTIONS TAKES A WHILE, SO IT'S NICE TO KEEP THIS COMMENTED OUT
# AND JUST LOAD THEM BELOW.
##step1. get the null distributions for counts of markers below a p theshold
##note that the p  values in the dataframe are not transformed, and these are used to get the window p values
# setwd('/Users/grovesdixon/Documents/lab_files/CpGoe_Project/trascriptome_only/data_analysis')
# reps = 10^5
# grab.size = 15
# cut = 0.05
# get.bootstrap = function(DF, id, grab.size, cut){
#   counts = c()
#   ps = c()
#   for (i in 1:reps){
#     sub = sample(DF$p, grab.size)
#     sub2 = sub[sub < cut]
#     new.count = length(sub2)
#     new.p = mean(sub)
#     counts = append(counts, new.count)
#     ps = append(ps, new.p)
#   }
#   results = data.frame(counts, ps)
#   write.table(results, paste('nullDist9-9-18', id, grab.size, cut, sep = "_"), row.names = F, quote = F)
#   return(results)
# }##Arugments: DF = the dataframe; id = a string to id the file from when exported
# ac.boot.df = get.bootstrap(ac.df, 'AC', grab.size, cut)
# ca.boot.df = get.bootstrap(ca.df, 'CA', grab.size, cut)
# con.boot.df = get.bootstrap(con.df, 'Control', grab.size, cut)
# hist(ca.boot.df$ps)
####################################################################################

#LOAD NULL DISTRIBUTIONS OF P VALUES
#These are permutations of random sets of 15 variants
#For each random set, the mean p value and the number of variants with
#unadjusted p values < 0.05 were counted. They provide the null
#distribution for clustering of moderately significant variants
grab.size = 15
cut = 0.05
ac.boot.df = read.table('nullDist9-9-18_AC_15_0.05', header = T)
ca.boot.df = read.table('nullDist9-9-18_CA_15_0.05', header = T)
con.boot.df = read.table('nullDist9-9-18_Control_15_0.05', header = T)

#We went with counts of variants below 0.05
#so set up vectors for those
ac.boot = ac.boot.df$counts
ca.boot = ca.boot.df$counts
con.boot = con.boot.df$counts


#SET UP THE FUNCTIONS FOR PERFORMING THE BOOTSTRAPPING

#FUNCTION get.boot.p
#A putatively selected variant is one with an unadjustd p value 0.05
#This function compares a count of putatively selected variants
#To the null distribution and returns the proportion of null entires 
#with equal more greater number of putatively selected variants
get.boot.p = function(null, obs){
  tot = length(null)
  more.extreme.count = length(null[null >= obs])
  p.value = more.extreme.count/tot
  return(p.value)
}


#FUNCTION setup.windows
#Divides the linkage groups into windows each with  number of markers equal to grab.size
setup.windows = function(grab.size, lg.tot){
	first.whole = floor(lg.tot/grab.size)##sets up the first windows up until the last
	left.sides = c()
	right.sides = c()
	left = 0
	for (i in 1:first.whole){
		right = grab.size*i
		left.sides = append(left.sides, left)
		right.sides = append(right.sides, right)
		left = right  ##now make the new left the previous right side, this way all adjacent windows will overlap
	}
	last.right = lg.tot
	last.left = lg.tot - (grab.size-1)
	left.sides = append(left.sides, last.left)
	right.sides = append(right.sides, last.right)
	window.sides = data.frame(left.sides, right.sides)
	return(window.sides)
}

#FUNCTION plot.boot.bars
#This combines the do.plot() function used above with the bootstrapping functions to plot
#a finalized figure that shows the combined p values and color coded bars below regions
#with significant bootstrapping signal
plot.boot.bars = function(DF, DF.ps, grab.size, cut, null.dist, MAIN, barYcoord, easy.cut, hard.cut){
  #ARGUMENTS: 
  # DF = the dataframe with the p values output for a particular cross (output from function do.plot()) 
  # DF.ps = data frame with the one-sided p values for each individual replicate
  # grab.size = the window size (in number of loci) used to build the null distribution beging used
  # cut = the cutoff p value cutoff for the null distribution
  # null.dist = the null distribution output from get.bootstrap()
  # MAIN = title
  # barYcoord = the vertical coordinate for placing the bootstrap bars
  # easy.cut = the unconservative p value cutoff for bootstrap bars
  # hard.cut = the conservative cutoff
  x = c()
  boot.ps = c()
  sig.counts = c()
  window.left.bounds = c()
  window.right.bounds = c()
  for (i in 1:14){
    lg.sub = DF[DF$LG == i,] ##pull out the linkage groups one at a time
	##set up the windows. This is more complicated than just seq(by = grab.size, because we want to include every marker)
	window.sides = setup.windows(grab.size, length(lg.sub[,1]))##sets up the windows as a dataframe of left and right sides
    for (win in seq(1:(length(window.sides[,1])))){
      left = window.sides[win,1]
      right = window.sides[win,2]
      sub = lg.sub[left:right,]
      win.left.cM = sub[1,'CM']##grabs the left side of this window in cM
      win.right.cM = sub[length(sub[,1]),'CM']##grabs the right side of the window
      window.left.bounds = append(window.left.bounds, win.left.cM)
      window.right.bounds = append(window.right.bounds, win.right.cM)
      #############################
      #get p value for counts of significant markers (turn this on and off)
      sig.sub = sub$p[sub$p < cut]
      obs.count = length(sig.sub)
      p = get.boot.p(null.dist, obs.count)
      #############################
      #############################
      #############################
      boot.ps = append(boot.ps, p)
      cM = median(sub$CM)
      x = append(x, cM)
    }
  }
  #This is the part where we trace the bars indicating regions with significant boostrapping signal
  Ycoords.easy = p.adjust(boot.ps, method = 'BH') ###set up the Y coordinates for the bars with the easier thrshold
  Ycoords.hard = p.adjust(boot.ps, method = 'BH') ##set up the Y coordiantes for hte barts with the harder threshold
  ##set up the Y coordinates for the 0.05 group of bars
  Ycoords.easy[Ycoords.easy < hard.cut] <- 100 ##don't mark the more significant ones here, so colors don't overlap
  Ycoords.easy[Ycoords.easy < easy.cut] <- barYcoord ##set the Y coordinates to the correct spot for the acutal significant ones by easy cutoff
  Ycoords.easy[Ycoords.easy != barYcoord] <- 100 ##set coordinate for all others outside so they won't be plotted
  ##set up the Y coordinates for the 0.05 adjusted group of bars
  Ycoords.hard[Ycoords.hard < hard.cut] <- barYcoord ##set coordiantes for significant hard cutoff bars
  Ycoords.hard[Ycoords.hard != barYcoord] <- -100 ##and get rid of those that arn't
  do.plot(DF.ps, MAIN, 100, 5, 100) ##do the plot
  par(lend = 3) ##set the option of how to end the line segments. 3 gives sharper ends.
  bottom = (barYcoord/3) ##this is to stack thin lines, there is a bottom and top for the range of Y coordiantes to stack the line segments in
  top = -barYcoord/3 ##the - sign makes it positive
  width = 2 ##set up the width of hte line to use
  for (Y in seq(bottom, top, by = 0.01)){
  	segments(window.left.bounds, Ycoords.easy + Y, window.right.bounds, Ycoords.easy+ Y, lwd = width, col = BAR.COLOR)
  } ###for the Y coordinates in the range between top and bottom, draw line segments that will make up the bars for hte easy cutoff
  for (Y in seq(bottom, top, by = 0.01)){
  segments(window.left.bounds, Ycoords.hard + Y, window.right.bounds, Ycoords.hard+ Y, lwd = width, col = BAR.HARD.COLOR)
  } ##draw the bars for the hard cutoff
  adj.boot.ps = p.adjust(boot.ps, method = 'BH')
  results = data.frame(x,Ycoords.easy, Ycoords.hard, boot.ps, adj.boot.ps, window.left.bounds, window.right.bounds)
  return(results)
}

#######################################################
#################### DO FINAL PLOT ####################
#######################################################
###SET UP SOME PLOTTING VARIABLES COLOR CODING SCHEME
quartz()
par(mfrow = c(3,1))
setYlim = c(-.2, 8)
BAR.COLOR = 'darkgreen'
BAR.HARD.COLOR = 'green'
HARD.CUT = 0.05
EASY.CUT = 0.1
#PLOT
ca.bars = plot.boot.bars(ca.df, ca.ps, grab.size, cut, ca.boot, 'Heat Selection in CA Cultures', -.25, EASY.CUT, HARD.CUT)
ac.bars = plot.boot.bars(ac.df, ac.ps, grab.size, cut, ac.boot, 'Heat Selection in AC Cultures', -.25, EASY.CUT, HARD.CUT)
con.bars = plot.boot.bars(con.df, con.ps, grab.size, cut, ac.boot, 'Control Comparisons', -.25, EASY.CUT, HARD.CUT)

################################################################################
########################## CALCULATE STRENGTH OF SELECTION #####################
################################################################################
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


#GET CONTROL FREQUENCY CHANGES
# ac.cont.freq = get.frequencies(dat, 'AC1c_tot', 'AC2c_tot', 'AC1c_refCount', 'AC2c_refCount', 'AC1c_altCount', 'AC2c_altCount')
# ca.cont.freq = get.frequencies(dat, 'CA1c_tot', 'CA2c_tot', 'CA1c_refCount', 'CA2c_refCount', 'CA1c_altCount', 'CA2c_altCount')
# head(ac.cont.freq)

###################################################################################
### LOOK AT CHANGE IN ALLELE FREQUENCY FOR PARTICULAR CHROMOSOMES ACROSS GENOME ###
###################################################################################

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
ac.delt.ps$fisher.p = ac.combined.p
ca.delt.ps$fisher.p = ca.combined.p


#ALSO UPDATE THE CONTROL COMPARISON
#here the delta p values were calculated as CAcontrol1 - CAcontrol2 and ACcontrol1 - ACcontrol2
con.delt.ps$variant = VARIANT
con.delt.ps$cM = CM
con.delt.ps$fisher.p = con.combined.p

#BUILD A CUMULATIVE PHASING DATAFRAME FROM THE PHASING DATA FOR EACH LINKAGE GROUP
for (i in lgs){
  phase.file = paste('g', i, '.loc.phasing.txt', sep = "")
  pdat = read.table(phase.file, header = T)
  pdat$lg = i
  endpoint = endpoints[i]
  pdat$cum.pos = pdat$position + endpoint
  pdat2 = merge(pdat, p.dat, by = 'variant')
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
# AC.control.dat$CA1dp = chrom.dp(CA.control.dat, 'C.c1', 'ref', 'C', 'ref.single', 'AC.Controls.delta.p')
# AC.control.dat$CA2dp = chrom.dp(CA.control.dat, 'C.c2', 'ref', 'C', 'ref.single', 'AC.Controls.delta.p')
#for the ca controls
# CA.control.dat$AC1dp = chrom.dp(CA.control.dat, 'A.c1', 'ref', 'A', 'ref.single', 'CA.Controls.delta.p')
# CA.control.dat$AC2dp = chrom.dp(CA.control.dat, 'A.c2', 'ref', 'A', 'ref.single', 'CA.Controls.delta.p')
CA.control.dat$CA1dp = chrom.dp(CA.control.dat, 'C.c1', 'ref', 'C', 'ref.single', 'CA.Controls.delta.p')
CA.control.dat$CA2dp = chrom.dp(CA.control.dat, 'C.c2', 'ref', 'C', 'ref.single', 'CA.Controls.delta.p')

#MAKE TWO SUBSET DATAFRAMES OF ONLY THE HETEROZYGOUS VARIANTS FOR EACH PARENT
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

#FUNCTION plot.chroms
#Plots the chromosome specific allele frequency changes for a given parent
plot.chroms = function(dat, chrom1.col, chrom2.col, YLIM, MAIN, XLAB, YLAB, COLORS, CUT){
  plot(dat[,chrom1.col]~dat$cM, type = "n", ylim = YLIM, main = MAIN, axes = F, xlab = XLAB, ylab = YLAB)
  lines(dat[,chrom1.col]~dat$cM, col = COLORS[1], lwd = 1)
  lines(dat[,chrom2.col]~dat$cM, col = COLORS[2], lwd = 1)
  abline(v = endpoints, col = 'grey', lwd = 0.5)
  sig = dat[p.adjust(dat$fisher.p, method = "BH") < CUT,]
  points(abs(sig[,chrom1.col])~sig$cM, pch = 19, col = 'red', cex = 0.5)
#   legend(1410, 0, c("Chromosome 1", "Chromosome 2"), fill = c('red', 'blue'))
  axis(1, at = endpoints, labels = F)
  axis(2, las = 1)
  mtext(lg.names, side = 1, at = mids, line = .75)
}

#PLOT FOR THE AC REPLICATES
YLIM = c(-.3, .3)
COLORS = c('cyan3', 'coral')
CUT = 0.05
quartz()
par(mfrow = c(2,1))
titleM = "AC (heat vs control) Maternal Chromosomes"
titleP = "AC (heat vs control) Paternal Chromosomes"
YLAB = expression(paste('P'['heat'], " - P"["control"], sep = ""))
plot.chroms(ac.a.poly, 'AC1dp', 'AC2dp', YLIM, titleM, "Linkage Group", YLAB, COLORS, CUT)
plot.chroms(ac.c.poly, 'CC1dp', 'CC2dp', YLIM, titleP, "Linkage Group", YLAB, COLORS, CUT)

#PLOT FOR THE CA REPLICATES
quartz()
par(mfrow = c(2,1))
titleM = "CA (heat vs control) Maternal Chromosomes"
titleP = "CA (heat vs control) Paternal Chromosomes"
plot.chroms(ca.c.poly, 'CC1dp', 'CC2dp', YLIM, titleM, "Linkage Group", YLAB, COLORS, CUT)
plot.chroms(ca.a.poly, 'AC1dp', 'AC2dp', YLIM, titleP, "Linkage Group", YLAB, COLORS, CUT)

#PLOT FOR THE CONTROLS
quartz()
par(mfrow = c(2,1))
titleM = "AC (control vs control) A Chromosomes"
titleP = "CA (control vs control) C Chromosomes"
plot.chroms(ac.control.poly, 'AC1dp', 'AC2dp', YLIM, titleM, "Linkage Group", YLAB, COLORS, CUT)
plot.chroms(ca.control.poly, 'CA1dp', 'CA2dp', YLIM, titleP, "Linkage Group", YLAB, COLORS, CUT)


#BOOTSTRAPPING WITH CHROMOSOME SPECIFIC ALLELE FREQUENCY CHANGES

bootstrap.chroms = function(iter, grab.size, dat, chrom.delta.ps){
  results = c()
  for (i in 1:iter){
    sample = mean(sample(abs(dat[,chrom.delta.ps]), grab.size), replace = F)
    results = append(results, sample)
  }
  return(results)
}

get.boot = function(lgs, dat, chrom.ds, grab.size, overlap.portion, boot.dat){
  means = c()
  positions = c()
  boot.pvalues = c()
  for (lg in lgs){
    lg.sub = dat[dat$lg == lg,]
    windows = seq(1, nrow(lg.sub), by = grab.size*overlap.portion)
    for (w in 1 : (length(windows) - 1)){
      left = windows[w]
      right = windows[(w+1)]
      win.sub = lg.sub[left:right,]
      mn = mean(abs(win.sub[,chrom.ds]))
      boot.p = length(boot.dat[boot.dat >= mn]) / length(boot.dat) + 1/length(boot.dat)
      means = append(means, mn)
      positions = append(positions, median(win.sub$cM))
      boot.pvalues = append(boot.pvalues, boot.p)
    }
  }
  results = data.frame(means, positions, boot.pvalues)
  colnames(results) = c(chrom.ds, 'cM', 'boot.pvalues')
  
  results$p = p.adjust(results$boot.pvalues, method = 'BH')
  sig = results[results$p < CUT,]
  sig$zero = 0
  points(sig$zero~sig$cM, pch = 15, col = 'red')
  sig = results[results$p < CUT2,]
  sig$zero = 0
  points(sig$zero~sig$cM, pch = 15, col = 'purple')
  return(results)
}


ITER = 10^5
ac.Ac1.boot = bootstrap.chroms(ITER, 30, ac.a.poly, 'AC1dp')
ac.Ac2.boot = bootstrap.chroms(ITER, 30, ac.a.poly, 'AC2dp')
ac.Cc1.boot = bootstrap.chroms(ITER, 30, ac.c.poly, 'CC1dp')
ac.Cc2.boot = bootstrap.chroms(ITER, 30, ac.c.poly, 'CC2dp')

ca.Ac1.boot = bootstrap.chroms(ITER, 30, ca.a.poly, 'AC1dp')
ca.Cc1.boot = bootstrap.chroms(ITER, 30, ca.c.poly, 'CC1dp')
con.Ac1.boot = bootstrap.chroms(ITER, 30, ac.control.poly, 'AC1dp')
con.Cc1.boot = bootstrap.chroms(ITER, 30, ca.control.poly, 'CA1dp')
hist(ac.Ac1.boot)


fraction = 1/10

#run 
plot.chroms(ac.c.poly, 'CC1dp', 'CC2dp', YLIM, 'AC paternal', "Linkage Group", YLAB, COLORS, CUT)
x = get.boot(lgs, ac.c.poly, 'CC1dp', 30, fraction, ac.Cc1.boot)
plot.chroms(ca.a.poly, 'AC1dp', 'AC2dp', YLIM, 'CA paternal', "Linkage Group", YLAB, COLORS, CUT)
x = get.boot(lgs, ca.a.poly, 'AC1dp', 30, fraction, ac.Ac1.boot)

#controls
plot.chroms(ac.control.poly, 'AC1dp', 'AC2dp', YLIM, 'AC control', "Linkage Group", YLAB, COLORS, CUT)
x = get.boot(lgs, ac.control.poly, 'AC1dp', 30, fraction, con.Ac1.boot)
plot.chroms(ca.control.poly, 'CA1dp', 'CA2dp', YLIM, 'CA control', "Linkage Group", YLAB, COLORS, CUT)
x = get.boot(lgs, ca.control.poly, 'CA1dp', 30, fraction, con.Cc1.boot)


x = get.boot(lgs, ac.c.poly, 'CC1dp', 30, 1/3, ac.Cc1.boot)
x = get.boot(lgs, ac.control.poly, 'AC1dp', 30, 1/3, con.Ac1.boot)
x = get.boot(lgs, ca.control.poly, 'CA1dp', 30, 1/3, con.Cc1.boot)
x = get.boot


x$p = p.adjust(x$boot.pvalues, method = 'BH')
plot(-log(x$p, 10)~x$cM)
abline(v = endpoints)

##################################################################
########## LOOKING AT WHICH GENES ARE BENEATH THE PEAKS ##########
##################################################################
setwd("./gene_data/")
#here the center point for the major peaks are identified as the variant with the highest -log(p)
#we then look for orthologs from the Acropora millepora transcriptome (Moya 2011) that align (Blast) near that locus


#UPLOAD ORTHOLOG (GENE) COORDINATES ON THE MAP
##This method works with the python script iso_to_LG.py and some blast results
genes = read.table('gene_lg_coordinates.txt', header = T, sep = "\t")
##This dataset is generated from the python script iso_to_LG.py and some blast results
##it includes all putative genes and was generated with this command: iso_to_LG.py -m linkage_map_Rinput.txt -s allIsogroups.txt -i amil_iso2gene.tab -o allGenes2LgOut.txt -a amil_Blast2_amil_genome_fold_c.out -tig2group amil_seq2iso.txt
##for more on gene annotation see Matz lab website: http://www.bio.utexas.edu/research/matz_lab/matzlab/Methods.html

#GET THE CENTER POINT FOR EACH GENE
genes$scaff.pos = apply(genes[3:4], 1, mean)

#FUNCTION get.scaff.pos()
#pulls scaffold positions from each variant name
#each variant name contains its scaffold coordinate (c486289_87218 is located on scaffold 486289 at position 87218)
get.scaff.pos = function(df, varcol, split.position){
  vec = c() ##place to store either bit of RAD tag
  for (i in df[,varcol]){
    dat = strsplit(i, '_')[[1]][split.position]
    vec = append(vec, dat)
  }
  return(vec)
}

#FUNCTION get.top()
#returns the most significant variant for a given linkage group
#if multiple points 
get.top = function(df, lg){
  df$adj.logp = -log(p.adjust(df$p, method = 'BH'), 10)
  sub = df[df$LG == lg,]
  sub = sub[sub$p == min(sub$p),]
  point = sub$CM
  print(paste("Top point =", point))
  return(sub)
}

#FUNCTION get.sub()
#pull only the significant variants for a given linkage group
get.sub = function(df, lg, cut){##subset a mapping dataframe for significant markers on a given lg
  df$adj.logp = -log(p.adjust(df$p, method = 'BH'), 10)
  sub = df[df$LG == lg,]
  sub = sub[sub$adj.logp > cut,]
  sub = sub[with(sub, order(VARIANT)),]
}

#SOME FUNCTIONS FOR PULLING GENES BASED ON DIFFERENT CRITERIA:

#FUNCTION pull.genes.based.on.lg()
##pulls genes close to a given map position in a given linkage group
pull.genes.based.on.lg = function(lg, pos, genes.df, win.size){
  lg.sub = genes.df[genes.df$lg == lg,]
  left = pos - (.5*win.size)
  right = pos + (.5*win.size)
  win.sub = lg.sub[lg.sub$cM > left,]
  win.sub = win.sub[win.sub$cM < right,]
  print(paste("Found", length(win.sub[,1]), "genes within the window"))
  return(win.sub)
}##pulls genes based on predicted linkage map position assigned as that of the closest RAD marker

#FUNCTION pull.genes.based.on.scaffold
#pulls isogroups that map to a specified genomic scaffold within a given range of a given point
#returns both annotated an unannotated isogroups
pull.genes.based.on.scaffold = function(scaff, genes.df, bp.left, bp.right){
  lg.sub = genes.df[genes.df$scaffold == scaff,]#subset the scaffold
  win.sub = lg.sub[lg.sub$scaff.pos > bp.left,]#subset for left and right bounds
  win.sub = win.sub[win.sub$scaff.pos < bp.right,]
  print(paste("Found", length(win.sub[,1]), "genes within the window"))
  return(win.sub)
}#pulls genes from a specified scaffold, within specified bp boundaries

#FUNCTION get.best.ranked()
#pulls the top X closest genes to a given position in a genomic scaffold
#only returns annotated genes
get.best.ranked = function(X, top.scaff, top.bp){
  g = genes[genes$scaffold == top.scaff,]
  g = na.omit(g)
  g$dist = abs(g$scaff.pos - top.bp)
  g = g[with(g, order(dist)),]
  best.ranked = g[1:X,]
  print(paste("Top", X, "Named Genes in proximity to top RAD marker"))
  return(best.ranked)
}


#PREP THE DATASETS
ca.df$scaff = get.scaff.pos(ca.df, 'VARIANT', 1)
ca.df$pos = get.scaff.pos(ca.df, 'VARIANT', 2)
ac.df$scaff = get.scaff.pos(ac.df, 'VARIANT', 1)
ac.df$pos = get.scaff.pos(ac.df, 'VARIANT', 2)

#########################################################
############ LOOK AT GENES NEAR PEAK AT LG 5 ############
#########################################################
#PULL THE FIGURE BACK UP
par(mfrow = c(1,1))
do.plot(ac.ps, 'Selection in AC Replicates', 100, 5, 100)
#ASSIGN CRITERIA FOR BEING NEAR THE PEAK
lg = 5
rank = 15
cM.win.size = 5
bp.win.size = 200000

#SUBSET THE DATASET
peak5 = get.sub(ac.df, lg, 2)
#PULL THE TOP LOCUS
top = get.top(ac.df, lg)
#ASSIGN THE COORDINATES FOR THE WINDOW WE ARE LOOKING IN
top.scaff = top$scaff
top.bp = as.numeric(top$pos)
top.cM = top$CM
abline(v = top.cM)
bp.left = top.bp - bp.win.size*.5
bp.right = top.bp + bp.win.size*.5

#SET UP THE GENES FOUND ON THIS SCAFFOLD
g = genes[genes$scaffold == top.scaff,]

#PULL THE GENES BASED ON LINKAGE GROUP
top.genes = pull.genes.based.on.lg(lg, top.cM, genes, cM.win.size)
print(top.genes)

#PULL THE GENES BASED ON SCAFFOLD
top.lg5.scaffold  = pull.genes.based.on.scaffold(top.scaff, genes, bp.left, bp.right)
print(top.lg5.scaffold )

#GRAB THE TOP GENES RANKED BY PROXIMITY TO CENTER OF PEAK
top15.lg5 = get.best.ranked(rank, top.scaff, top.bp)


##########################################################
############ LOOK AT GENES NEAR PEAK AT LG 10 ############
##########################################################
#LOOK AT THE PEAK
do.plot(ca.ps, 'Selection in CA Replicates', 100, 5, 100)
#SUBSET FOR THE PEAK
peak10 = get.sub(ca.df, 10, 2)
#SET PROXIMITY CIRTERIA
lg = 10
cut = 2
rank = 15
cM.win.size = 3
bp.win.size = 200000
bp.left = top.bp - bp.win.size*.5
bp.right = top.bp + bp.win.size*.5

#PULL THE TOP MARKER AND ITS COORDINATES
top = get.top(ca.df, lg)
top.bp = as.numeric(top$pos)
top.cM = top$CM
abline(v = top.cM)
top.scaff = top$scaff

#GET THE LISTS OF GENES
g = genes[genes$scaffold == top.scaff,]##the genes that align to the top scaffold
top = get.top(ca.df, lg)
top.bp = as.numeric(top$pos)
bp.left = top.bp - bp.win.size*.5
bp.right = top.bp + bp.win.size*.5
top.lg10.lg = pull.genes.based.on.lg(lg, top$CM, genes, cM.win.size) ##get the genes based on predicted locations
top.lg10.scaffold = pull.genes.based.on.scaffold(top.scaff, genes, bp.left, bp.right)##get the genes based on scaffold
top15.lg10 = get.best.ranked(rank, top.scaff, top.bp) ## get the top ranked genes in proximity to the most significant RAD locus (only includes annoated)

print(top.lg10.lg)
print(top.lg10.scaffold)
print(top15.lg10)

#WRITE OUT THE GENES FROM LG 1- PEAK
write.table(top.lg10.scaffold, "./genes_under_lg10_peak.txt", quote = F, row.names = F)



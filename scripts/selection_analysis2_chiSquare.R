#RUN CHI SQUARE ANALYSIS ON ALLELE FREQUENCY DIFFERENCES AND PLOT MANHATTAN PLOTS

#LOAD THE DATA OUTPUT FROM selection_analysis_setup.R
path.2.working.dir = "~/git_Repositories/linkage_mapping_and_selection/data_files/"
setwd(path.2.working.dir)
load('LinkageMappingSetup.R')

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

#FUNCTION plot.scan()
#Adjusts a set of p values for FDR and plots them based on map locations and adds bars below to indicate boostrap significances
#ARGUMENTS: 
# DF = the dataframe
# MAIN = the title you want for the plot
# LEFT = vector of the left sides of bootstrap bars
# Y = the Y coordinate for where to plot the bootstrap bars
# LEN = vector of lengths for the bootstrap bars
plot.scan = function(DF, MAIN, LEFT, Y, LEN){##NOTE THIS FUNCTION IS REPEATED BELOW
  CEX = 1.2
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
  plot(-log(p, 10)~CM, data = DF, main = NULL, pch = 19, col = colors, cex = point.size, axes = F, xlab = "Linkage Group", ylab = expression('-log'['10']*'(p)'), ylim=setYlim)
  axis(1, at = endpoints, labels = F)
  axis(2, labels = NULL, las = 1, cex.axis = CEX)
  mtext(lg.names, side = 1, at = mids, line = .75, cex = CEX)
  par(lend = 2)
  if (PLOT.SCALE == TRUE){
    segments(LEFT, Y, (LEFT + LEN), Y, lwd = 2)
    TEXT = paste(LEN, "cM")
    labPos = LEFT + .5 * LEN
    text(labPos, y = (Y + .4), labels = TEXT)
  }
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

## NOW WE CAN PLOT THE MANHATTAN PLOTS FOR EACH REPLICATE
####### DO AC ONE-TAILED PLOT ##############
ac1.one.side = get.one.side.ps(dat, "AC1c_tot", "AC1h_tot", "AC1c_refCount", "AC1h_refCount", "AC1c_altCount", "AC1h_altCount", ac.alternative)
ac2.one.side = get.one.side.ps(dat, "AC2c_tot", "AC2h_tot", "AC2c_refCount", "AC2h_refCount", "AC2c_altCount", "AC2h_altCount", ac.alternative)
ac.ps = data.frame(ac1.one.side, ac2.one.side)
head(ac.ps)
par(mfrow = c(1,1))
PLOT.SCALE = TRUE
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
  abline(h = 0.25)
}
quartz()
LOESS.SPAN = .05
par(mfrow = c(4,1))
plot.abs(ac.delt.ps, 'ac1.ps.delt.p', 'AC1 changes in allele frequencies')
plot.abs(ac.delt.ps, 'ac2.ps.delt.p', 'AC2 changes in allele frequencies')
plot.abs(ca.delt.ps, 'ca1.ps.delt.p', 'CA1 changes in allele frequencies')
plot.abs(ca.delt.ps, 'ca2.ps.delt.p', 'CA2 changes in allele frequencies')
plot.abs(con.delt.ps, 'AC.Controls.delta.p', 'Control Changes in allele frequencies')
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
#with significant bootstrapping signal. Also outputs the bootstrap coordiantes for later use.
plot.boot.bars = function(DF, DF.ps, grab.size, cut, null.dist, MAIN, bar.Ycoord, bar.width, easy.cut, hard.cut){
  #ARGUMENTS: 
  # DF = the dataframe with the p values output for a particular cross (output from function do.plot()) 
  # DF.ps = data frame with the one-sided p values for each individual replicate
  # grab.size = the window size (in number of loci) used to build the null distribution beging used
  # cut = the cutoff p value cutoff for the null distribution
  # null.dist = the null distribution output from get.bootstrap()
  # MAIN = title
  # barYcoord = the vertical coordinate for placing the bootstrap bars
  # bar.width = the vertical width of the bootstrap bar
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
      #get p value for counts of significant markers (turn this on and off)
      sig.sub = sub$p[sub$p < cut]
      obs.count = length(sig.sub)
      p = get.boot.p(null.dist, obs.count)
      boot.ps = append(boot.ps, p)
      cM = median(sub$CM)
      x = append(x, cM)
    }
  }
  result1 = data.frame(x, boot.ps, window.left.bounds, window.right.bounds)
  result1$adj.ps = p.adjust(result1$boot.ps, method = "BH")
  easy = result1[result1$adj.ps > hard.cut,]
  easy = easy[easy$adj.ps < easy.cut,]
  hard = result1[result1$adj.ps < hard.cut,]
  print(head(easy))
  print(nrow(easy))
  do.plot(DF.ps, MAIN, 100, 5, 100)
  par(xpd=TRUE)
  top = bar.Ycoord + 0.5*bar.width
  bottom = bar.Ycoord - 0.5*bar.width
  print(seq(bottom, top, by = 0.01))
  for (y in seq(bottom, top, by = 0.01)){
    if (nrow(easy) >= 1){
      for (i in 1:nrow(easy)){
        segments(easy$window.left.bounds[i], y, easy$window.right.bounds[i], y, col = EASY.BAR.COLOR, lwd = 2)
      }
    }
    if (nrow(hard) >= 1){
      for (i in 1:nrow(hard)){
        segments(hard$window.left.bounds[i], y, hard$window.right.bounds[i], y, col = HARD.BAR.COLOR, lwd = 2)
      }
    }
  }
  par(xpd=FALSE)
  if (nrow(hard) > 0){
    if (nrow(easy) > 0){
      print("easy:")
      print(easy)
      print("hard:")
      print(hard)
      e.lefts = easy$window.left.bounds
      e.rights = easy$window.right.bounds
      h.lefts = hard$window.left.bounds
      h.rights = hard$window.right.bounds
      output.e = data.frame(e.lefts, e.rights)
      output.e$cut = "easy"
      output.h = data.frame(h.lefts, h.rights)
      output.h$cut = "hard"
      COLNAMES = c('lefts', 'rights', 'cut')
      colnames(output.e) = COLNAMES
      colnames(output.h) = COLNAMES
      output = rbind(output.e, output.h)
      print(output)
      return(output)
    }
  }
}

#######################################################
#################### DO FINAL PLOT ####################
#######################################################
###SET UP SOME PLOTTING VARIABLES COLOR CODING SCHEME
quartz()
par(mfrow = c(1,1))
EASY.BAR.COLOR = 'darkgreen'
HARD.BAR.COLOR = 'green'
HARD.CUT = 0.05
EASY.CUT = 0.1
PLOT.SCALE = TRUE
#PLOT
setYlim = c(-.75, 8)
point.size = .75
par(mar=c(2,2,1,2)+0.1)
ca.bars = plot.boot.bars(ca.df, ca.ps, grab.size, cut, ca.boot, 'Heat Selection in CA Cultures', .825*setYlim[1], -.75*setYlim[1], EASY.CUT, HARD.CUT)
setYlim = c(-.6, 6)
PLOT.SCALE = FALSE
par(mar=c(2,2,1,2)+0.1)
ac.bars = plot.boot.bars(ac.df, ac.ps, grab.size, cut, ac.boot, 'Heat Selection in AC Cultures', .8*setYlim[1], -.75*setYlim[1], EASY.CUT, HARD.CUT) #5.25
con.bars = plot.boot.bars(con.df, con.ps, grab.size, cut, ac.boot, 'Control Comparisons', .8*setYlim[1], -.75*setYlim[1], EASY.CUT, HARD.CUT)


################################################################################################
#################### OUTPUT THE BOOSTRAP BAR COORDINATES AND SAVE WORKSPACE ####################
################################################################################################
write.table(ca.bars, 'ca_bars.txt')
write.table(ac.bars, 'ac_bars.txt')

#OUTPUT THE WORKSPACE
save.image("selection_analysis2_output.R")



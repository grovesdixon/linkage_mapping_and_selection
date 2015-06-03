#FIND GENES POSITIONED NEAR SELECTION PEAKS

#LOAD THE DATA OUTPUT FROM selection_analysis2_chiSquare.R
path.2.working.dir = "~/git_Repositories/linkage_mapping_and_selection/data_files/"
setwd(path.2.working.dir)
load('selection_analysis2_output.R')

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

print('Genes Near Selection Peak on Linkage Group 10:')
print(top.lg10.scaffold)

#WRITE OUT THE GENES FROM LG 1- PEAK
write.table(top.lg10.scaffold, "./genes_under_lg10_peak.txt", quote = F, row.names = F)

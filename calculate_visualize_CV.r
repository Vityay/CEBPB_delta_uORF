 
#####################################################################################
# Simple CV calculation
# by: Tristan de Jong
# Last edit: 17/01/2018
#####################################################################################
library("gProfileR")
library("EDASeq")                                        

# 0. CV calculation
getCV = function(makeCVData){
  SD1 = apply(makeCVData, 1, sd, na.rm=TRUE) 
  mean1 = rowMeans(makeCVData, na.rm = TRUE)
  CV1 = SD1/mean1
  return(CV1)
}

myColors = c("red","pink","blue","lightblue")


pdf("results/analysis.pdf")

##########################################################################
#1. Read the data
##########################################################################
rawData = read.delim("data/20170306_RNASeq_Gertrud_rawcounts_0FPM.txt", header =T,
                     sep= "\t", row.names = "Gene")

phenotypes = read.delim("data/Phenotype.csv", header =T, sep = ",", row.names = "Library")
phenotypes$AgeMut = paste(phenotypes$Age,phenotypes$Mutant,sep = "_")


##########################################################################
#2. Filter to a minimum of 1 FPM and remove 0 counts
##########################################################################
Filter0 = apply(rawData,1,min) >= 1                        # All data must be non-zero
rawData = rawData[Filter0,]                                # Apply the non-zero filter

pre_fpm <- sweep(rawData,2, colSums(rawData),'/')          # Divide by all collumns
FPM_normalized <- sweep(pre_fpm, 2, 1e6, '*')            # Multiply by 1.000.000
# FPM_filter =   rowMeans(FPM_normalized[,phenotypes$AgeMut == "young_wildtype" ]) >= 1 |
#                rowMeans(FPM_normalized[,phenotypes$AgeMut == "old_wildtype"   ]) >= 1 |
#                rowMeans(FPM_normalized[,phenotypes$AgeMut == "young_knock-in" ]) >= 1 |
#                rowMeans(FPM_normalized[,phenotypes$AgeMut == "old_knock-in"   ]) >= 1
FPM_filter =  rowMeans(FPM_normalized) >= 1
filteredData = rawData[FPM_filter,]    

# UQ normalize
UQ_normalized = betweenLaneNormalization(x = as.matrix(filteredData), which=c("upper"), offset=FALSE, round=FALSE)
boxplot(UQ_normalized, outline = FALSE, las = 2, cex.axis = 0.5, col = myColors[as.factor(phenotypes$AgeMut)])


##########################################################################
# 3. UQ normalized
##########################################################################
YW = UQ_normalized[,phenotypes$AgeMut == "young_wildtype"]
OW = UQ_normalized[,phenotypes$AgeMut == "old_wildtype"]
YK = UQ_normalized[,phenotypes$AgeMut == "young_knock-in"]
OK = UQ_normalized[,phenotypes$AgeMut == "old_knock-in"]

##########################################################################
#4 Get CVs
##########################################################################
YoungWT = getCV(YW)
YoungKI = getCV(YK)
OldWT = getCV(OW)
OldKI = getCV(OK)

dat = cbind(YoungWT, YoungKI,OldWT,OldKI)

dat = DataFrame(dat)


##########################################################################
# 5. CV vs CV
##########################################################################
par(mfrow=c(2,2))

myColors = c("#1B7B34","#1FB58F","#EAB126","grey","darkgrey","#F24C4E","forestgreen","lightblue","red","magenta")
myColors2 = c("#1B7B34","#1B7B3480",
              "#1FB58F","#1FB58F80",
              "#EAB126",  "#EAB12680",
              "grey", "#BEBEBE80",
              "darkgrey","#A9A9A980",
              "#F24C4E","#F24C4E80",
              "forestgreen", "#228B2280",
              "lightblue", "#0000FF80",
              "red","#FF000080",
              "magenta","#FF00FF80")
black = adjustcolor( "black", alpha.f = 0.1)


# A. Young
plot(dat$YoungWT, dat$YoungKI, main="CV YoungWT vs YoungKI",
     col=black, pch=16, cex = 0.4,
     ylim=c(0,2), xlim=c(0,2), 
     xlab='CV, YoungWT', ylab='CV, YoungKI')
#abline(a=0,b=1)
abline(lm(dat$YoungKI~dat$YoungWT+0),col='red', lty=2)

top = table(dat$YoungWT>dat$YoungKI)[1]
bot = table(dat$YoungWT>dat$YoungKI)[2]
#legend(x = "topleft",    legend = top, bty = "n")
#legend(x = "bottomright",legend = bot, bty = "n")


# B. Old
plot(dat$OldWT, dat$OldKI, main="CV OldWT vs OldKI",
     col=black, pch=16, cex = 0.4,
     ylim=c(0,2), xlim=c(0,2), 
      xlab='CV, OldWT', ylab='CV, OldKI')
#abline(a=0,b=1)
abline(lm(dat$OldKI~dat$OldWT+0),col='red', lty=2)

top = table(dat$OldWT>dat$OldKI)[1]
bot = table(dat$OldWT>dat$OldKI)[2]
#legend(x = "topleft",    legend = top, bty = "n")
#legend(x = "bottomright",legend = bot, bty = "n")

# C. WT
plot(dat$YoungWT, dat$OldWT, main="CV YoungWT vs OldWT",
     col=black, pch=16, cex = 0.4,
     ylim=c(0,2), xlim=c(0,2), 
     xlab='CV, YoungWT', ylab='CV, OldWT')
#abline(a=0,b=1)
abline(lm(dat$OldWT~dat$YoungWT+0),col='red', lty=2)

top = table(dat$YoungWT>dat$OldWT)[1]
bot = table(dat$YoungWT>dat$OldWT)[2]
#legend(x = "topleft",    legend = top, bty = "n")
#legend(x = "bottomright",legend = bot, bty = "n")

# D. KI
plot(dat$YoungKI, dat$OldKI, main="CV YoungKI vs OldKI",
     col=black, pch=16, cex = 0.4,
     ylim=c(0,2), xlim=c(0,2), 
     xlab='CV, YoungKI', ylab='CV, OldKI')
#abline(a=0,b=1)
abline(lm(dat$OldKI~dat$YoungKI+0),col='red', lty=2)

top = table(dat$YoungKI>dat$OldKI)[1]
bot = table(dat$YoungKI>dat$OldKI)[2]
#legend(x = "topleft",    legend = top, bty = "n")
#legend(x = "bottomright",legend = bot, bty = "n")


##########################################################################
#6. Boxplots
##########################################################################
par(mfrow=c(1,1))

boxplot(dat$YoungWT,dat$YoungKI,dat$OldWT, dat$OldKI,col = myColors,
        main = "Boxplot of BCV's",
        names = c("Young WT","Young KI","Old WT","Old KI"), las = 2)


##########################################################################
#7. Find DE
##########################################################################

## Effect of age on young
YoungWT_YoungKI = rownames(dat[dat$YoungWT  > dat$YoungKI *2,]) # Increase in noise
YoungKI_YoungWT = rownames(dat[dat$YoungWT  < dat$YoungKI /2,]) # Decrease in noise
Young_even      = rownames(dat[dat$YoungWT  < dat$YoungKI *2 & dat$YoungWT  > dat$YoungKI /2,]) # No change

# Effect of Knock in on old
OldWT_OldKI = rownames(dat[dat$OldWT  > dat$OldKI *2,]) # Increase in noise
write.csv(OldWT_OldKI,"OldWT_OldKI.csv")

OldKI_OldWT = rownames(dat[dat$OldWT  < dat$OldKI /2,]) # Decrease in noise
write.csv(OldKI_OldWT,"OldKI_OldWT.csv")

Old_even    = rownames(dat[dat$OldWT  < dat$OldKI *2 & dat$OldWT  > dat$OldKI /2,]) # No change

# Effect of age on noise
YoungWT_OldWT = rownames(dat[dat$YoungWT > dat$OldWT *2,]) # Increase in noise
OldWT_YoungWT = rownames(dat[dat$YoungWT < dat$OldWT /2,]) # Decrease in noise
WT_even       = rownames(dat[dat$YoungWT < dat$OldWT *2 & dat$YoungWT > dat$OldWT /2,]) # No change

# Effect of age on KI
YoungKI_OldKI = rownames(dat[dat$YoungKI > dat$OldKI *2,]) # Increase in noise
OldKI_YoungKI = rownames(dat[dat$YoungKI < dat$OldKI /2,]) # Decrease in noise
KI_even       = rownames(dat[dat$YoungKI < dat$OldKI *2 & dat$YoungKI > dat$OldKI /2,]) # No change

dev.off()

##########################################################################
#8. GO-terms
##########################################################################
pdf("results/GO_Terms.pdf")
par(mar=c(5.1,18.1,4.1,2.1), xpd = TRUE,mfrow=c(1,1)) #  bottom, left, top and right

######################################################################################################################################################
# # # Gprofile background
# background = rownames(dat)
# localGprofile   =   gprofiler(background, organism = "mmusculus", ordered_query = F,
#                               significant = F, exclude_iea = F, underrep = F, evcodes = F,
#                               region_query = F, max_p_value = 1, min_set_size = 0, max_set_size = 0,
#                               min_isect_size = 0, correction_method = "analytical",
#                               hier_filtering = "none", domain_size = "annotated", custom_bg = "",
#                               numeric_ns = "", png_fn = NULL, include_graph = F, src_filter = c("GO:BP", "GO:MF", "GO:CC"))
# 
# # src_filter = c("GO:BP", "GO:MF", "GO:CC")
# # src_filter = NULL
# #localGprofile = localGprofile[,1:13]
# write.table(localGprofile, file =  "data/GObackground2.csv", sep = "\t", row.names = FALSE)
######################################################################################################################################################

query = YoungWT_OldWT

goProfileThis = function(query, outname){
    localGprofile   =   gprofiler(query, organism = "mmusculus", ordered_query = F,
                                  significant = T, exclude_iea = F, underrep = F, evcodes = F,
                                  region_query = F, max_p_value = 1, min_set_size = 0, max_set_size = 0,
                                  min_isect_size = 0, correction_method = "analytical",
                                  hier_filtering = "none", #("none","moderate")
                                  domain_size = "annotated", custom_bg = "",
                                  numeric_ns = "", png_fn = NULL, include_graph = F, src_filter = c("GO:BP", "GO:MF", "GO:CC"))     # src_filter = c("GO:BP", "GO:MF", "GO:CC")
    #(GO:BP, GO:MF, GO:CC to select a particular GO branch), KEGG, REAC, TF, MI, CORUM, HP, HPA, OMIM. 
    # src_filter = NULL
    localGprofile = localGprofile[,1:13]
    rownames(localGprofile) =  localGprofile$term.name
    
    outFolder = "results/gprofiled/"
    outFile = paste(outname,".csv",sep = "")
    write.csv(localGprofile, file =  paste(outFolder,outFile, sep = ""), row.names = FALSE)
    return(localGprofile)
    }

plotGO = function(inSet, plotLabel){
    cutoff = 0.8
    a = reference[rownames(inSet),6]/reference[rownames(inSet),5] # GO terms for the reference
    b = inSet[,6]/inSet[,5]                                       # GO terms for my querie     
    
    b/a
    Interest = rownames(inSet[which(b/a > cutoff),])
    localCol = ifelse (b/a > cutoff, "red","black")
    range = c(min(c(a,b)), max(c(a,b)))
    
    if (length(Interest) >0){
      print(paste( plotLabel,length(Interest)))
      myBar = barplot(( inSet[Interest,6]/inSet[Interest,5]) /
                        (reference[Interest,6]/reference[Interest,5]), las = 2 ,
                      horiz = TRUE,
                      col = c("grey",myColors))
      axis(2, at=myBar[,1], 
           labels=inSet[Interest,12], 
           las = 2,
           cex.axis = ifelse(length(Interest)>10,0.4,0.8))
      
      text(inSet[Interest,6], y = myBar[,1], x = 0.4, cex = 0.5)
      
      title(main = paste(plotLabel,"\n","overrepresented"),
            cex.main = 0.75)
    }
    abline(v = 1, lty = 2)
    abline(v = 1.5, lty = 2)
    return(Interest)
}

# Load the background
reference = read.csv("data/GObackground.csv", sep = ",")
rownames(reference) = reference$term.name

# 1. Consequences of ageing without mutation.
GO  = goProfileThis(YoungWT_OldWT, "YoungWT_OldWT")
plotGO(GO, "YoungWT_OldWT")

GO  = goProfileThis(OldWT_YoungWT, "OldWT_YoungWT")
plotGO(GO, "YoungWT_OldWT")

# 2. Consequences of ageing with mutation
GO  = goProfileThis(YoungKI_OldKI , "YoungKI_OldKI")
plotGO(GO, "YoungKI_OldKI")

GO  = goProfileThis(OldKI_YoungKI, "OldKI_YoungKI")
plotGO(GO, "OldKI_YoungKI")

# 3. Consequences of ageing regardless of Mutation
twoGo = intersect(YoungWT_OldWT,YoungKI_OldKI)
GO  = goProfileThis(twoGo, "Both Young _ Both Old")
plotGO(GO, "Both Young _ Both Old")

twoGo = intersect(OldWT_YoungWT,OldKI_YoungKI)
GO  = goProfileThis(twoGo, "Both Old _ Both Young")
plotGO(GO, "Both Old _ Both Young")

# 4. Consequences mutation on Young
GO  = goProfileThis(YoungWT_YoungKI , "YoungWT_YoungKI")
plotGO(GO, "YoungWT_YoungKI")

GO  = goProfileThis(YoungKI_YoungWT , "YoungKI_YoungWT")
plotGO(GO, "YoungKI_YoungWT")

# 5. Consequences mutation on Old
GO  = goProfileThis(OldWT_OldKI, "OldWT_OldKI")
plotGO(GO, "OldWT_OldKI")

GO  = goProfileThis(OldKI_OldWT, "OldKI_OldWT")
plotGO(GO, "OldKI_OldWT")

# 6. The effect of KI regardless of age
twoGo = intersect(OldWT_OldKI,YoungWT_YoungKI)
GO  = goProfileThis(twoGo, "Both WT _ Both KI")
plotGO(GO,  "Both WT _ Both KI")

twoGo = intersect(OldKI_OldWT,YoungKI_YoungWT)
GO  = goProfileThis(twoGo, "Both KI _ Both WT")
plotGO(GO,  "Both KI _ Both WT")

dev.off()


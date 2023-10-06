
install.packages("CMplot")
library("CMplot")
library(data.table)
library(dplyr)


data(pig60K)   #calculated p-values by MLM
head(pig60K)
typeof(pig60K)


all<-fread("COG_EDU_SCZ_PLEEIO_Analysis_Results.txt.gz")
all<-all[,c("CHROM","POS","SNP","P")]
colnames(all)<-c("CHR","BP","SNP","P")
all2<-all[,c("CHR","BP","SNP","P")]
all2$PHE<-"ALL"
colnames(all)<-c("CHR","BP","SNP","ALL") 

dis<-fread("FUMA_INPUT_withoutMHC_Discordant.txt")
dis<-dis[, c("CHR","BP","SNP","P") ]
dis2<-dis[, c("CHR","BP","SNP","P") ]
dis2$PHE<-"DIS"
colnames(dis)<-c("CHR","BP","SNP","DIS") 

con<-fread("FUMA_INPUT_withoutMHC_Concordant.txt")
con<-con[, c("CHR","BP","SNP","P") ]
con2<-con[, c("CHR","BP","SNP","P") ]
con2$PHE<-"CON"
colnames(con)<-c("CHR","BP","SNP","CON") 

both<-fread("FUMA_INPUT_withoutMHC_Both.txt.gz")
both<-both[, c("CHROM","POS","SNP","PLEIO_META_PLEIO_P") ]
both2<-both[, c("CHROM","POS","SNP","PLEIO_META_PLEIO_P") ]
colnames(both2)<-c("CHR","BP","SNP","P") 
both2$PHE<-"BOTH"
colnames(both)<-c("CHR","BP","SNP","BOTH") 

# Perform outer merge on 'CHR', 'BP', and 'SNP'
merged_data <- merge(all, dis, by = c("CHR", "BP", "SNP"), all = TRUE)
merged_data <- merge(merged_data, con, by = c("CHR", "BP", "SNP"), all = TRUE)
merged_data <- merge(merged_data, both, by = c("CHR", "BP", "SNP"), all = TRUE)

merged_data<-merged_data[,c( "SNP","CHR","BP","ALL","DIS","CON","BOTH") ]
colnames(merged_data)<-c( "SNP","chr","pos","ALL","DIS","CON","BOTH")



merged_data<-as.data.frame(merged_data)
#merged_data <- merged_data %>%mutate_at(vars(ALL, DIS, CON), ~coalesce(., 0.999))
merged_data <- merged_data %>%mutate(ALL = ifelse(ALL == 0, 0.99, ALL))
merged_data <- merged_data %>%mutate(DIS = ifelse(DIS == 0, 0.99, DIS))
merged_data <- merged_data %>%mutate(CON = ifelse(CON == 0, 0.99, CON))
merged_data <- merged_data %>%mutate(BOTH = ifelse(BOTH == 0, 0.99, BOTH))

SNPs <-  merged_data[
	merged_data$ALL < 1e-8 |
	merged_data$DIS < 1e-8 |
    merged_data$CON < 1e-8 |
	merged_data$BOTH < 1e-8, 1]

CMplot(merged_data,type="p",plot.type="m",LOG10=TRUE,highlight.type="p",highlight=SNPs,
        threshold=1e-8,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
        signal.cex=0.5, signal.col="red", highlight.col="grey",highlight.cex=0.5,
        file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

merged_data2<-merged_data[,c("SNP","chr","pos","DIS","CON")]
CMplot(merged_data2, plot.type="m",multraits=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
        threshold.lwd=c(1,1), amplify=TRUE,bin.size=1e6,
        signal.cex=0.5, file="jpg",file.name="",dpi=300,file.output=TRUE,verbose=TRUE,
        highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1.4, 
        trait.legend.ncol=1, trait.legend.pos="left")  #signal.col=c("red","green"),chr.den.col=c("darkgreen", "yellow", "red"),threshold.col=c("black","grey")



######AJHG CONCORDANT Discordant

dis<-fread("ASSET_DISCORDANT.gz")
dis<-dis[,c("CHR","BP","SNP","P")]
colnames(dis)<-c("CHR","BP","SNP","AJHG_DIS") 

con<-fread("ASSET_CONCORDANT.gz")
con<-con[, c("CHR","BP","SNP","P") ]
colnames(con)<-c("CHR","BP","SNP","AJHG_CON") 

all<-fread("ASSET_ALL.txt")
all<-all[, c("CHR","BP","SNP","P") ]
colnames(all)<-c("CHR","BP","SNP","AJHG_ALL") 


# Perform outer merge on 'CHR', 'BP', and 'SNP'
merged_data2 <- merge(con, dis, by = c("CHR", "BP", "SNP"), all = TRUE)
merged_data2 <- merge(merged_data2, all, by = c("CHR", "BP", "SNP"), all = TRUE)


merged_data2<-merged_data2[,c( "SNP","CHR","BP","AJHG_CON","AJHG_DIS","AJHG_ALL") ]
colnames(merged_data2)<-c( "SNP","chr","pos","AJHG_CON","AJHG_DIS","AJHG_ALL")


merged_data2<-as.data.frame(merged_data2)
merged_data2 <- merged_data2 %>%mutate_at(vars(AJHG_DIS, AJHG_CON,AJHG_ALL), ~coalesce(., 0.999))
merged_data2 <- merged_data2 %>%mutate(AJHG_DIS = ifelse(AJHG_DIS == 0, 0.99, AJHG_DIS))
merged_data2 <- merged_data2 %>%mutate(AJHG_CON = ifelse(AJHG_CON == 0, 0.99, AJHG_CON))
merged_data2 <- merged_data2 %>%mutate(AJHG_CON = ifelse(AJHG_ALL == 0, 0.99, AJHG_ALL))


SNPs <-  merged_data2[
	merged_data2$AJHG_DIS < 1e-8 |
	merged_data2$AJHG_CON < 1e-8 |
	merged_data2$AJHG_ALL < 1e-8, 1]


CMplot(merged_data2,type="p",plot.type="m",LOG10=TRUE,highlight.type="p",highlight=SNPs,
        threshold=1e-8,threshold.col="black",threshold.lty=1,col=c("grey60","#4197d8"),
        signal.cex=0.5, signal.col="red", highlight.col="grey",highlight.cex=0.5,
        file="jpg",dpi=300,file.output=TRUE,verbose=TRUE,multracks=TRUE)

CMplot(merged_data2, plot.type="m",multraits=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
        threshold.lwd=c(1,1), amplify=TRUE,bin.size=1e6,
        signal.cex=0.5, file="jpg",file.name="",dpi=300,file.output=TRUE,verbose=TRUE,
        highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1.4, 
        trait.legend.ncol=1, trait.legend.pos="left")  #signal.col=c("red","green"),chr.den.col=c("darkgreen", "yellow", "red"),threshold.col=c("black","grey")







#### https://github.com/anastasia-lucas/hudson

library(hudson)
result <- rbind(all2, dis2, con2, both2)
colnames(result) <- c("CHR", "POS", "SNP", "pvalue", "PHE")
result <- result %>%mutate(pvalue = ifelse(pvalue == 0, 0.99, pvalue))
result<-result[,c( "PHE","SNP","CHR", "POS","pvalue")]

upper<-result[result$PHE=="ALL",]
lower <- result[(result$PHE != "ALL"), ]
upper<-as.data.frame(upper)
lower<-as.data.frame(lower)

#groupcolors=c(ALL="gray",DIS="green",CON="blue",BOTH="red")

phemirror(top=upper, bottom = lower, toptitle = "PheWAS Example: Data 1",opacity=0.1,groupcolors=c(ALL="dark gray",CON="red",DIS="yellow",BOTH="blue"), 
bottomtitle = "PheWAS Example: Data 2")



all2<-all[,c("SNP","CHR","BP","ALL")]
colnames(all2)<-c( "SNP","CHR","POS","pvalue")

dis2<-dis[,c("SNP","CHR","BP","DIS")]
colnames(dis2)<-c( "SNP","CHR","POS","pvalue")

con2<-con[,c("SNP","CHR","BP","CON")]
colnames(con2)<-c( "SNP","CHR","POS","pvalue")


both2<-con[,c("SNP","CHR","BP","CON")]
colnames(both2)<-c( "SNP","CHR","POS","pvalue")




gmirror(top=dis2, bottom=con2, tline=0.00000005, bline=0.00000005, 
toptitle="Discordant", bottomtitle = "Concordant", 
highlight_p = c(0.00000005,0.00000005), highlighter="#d70909")


all2$PHE<-"ALL"
dis2$PHE<-"DIS"
con2$PHE<-"CON"

dis_con <- rbind(dis2, con2)


phemirror(top=all2, bottom = dis_con, toptitle = "PheWAS Example: Data 1", 
bottomtitle = "PheWAS Example: Data 2")

##AJH concordant Discordant


##############---------------------------gwaslab-------------------------------
import sys 
sys.path.insert(0,"/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Inputs/FUMA/")
import gwaslab as gl
mysumstats = gl.plot_miami(path1="FUMA_INPUT_withoutMHC_Discordant.tsv" ,
                           path2="FUMA_INPUT_withoutMHC_Concordant.tsv",
                           cols1=["CHR","BP","P"],
                           cols2=["CHR","BP","P"],
                           titles=["bmi male","bmi female"],
                           titles_pad=[0.15,0.0],
                           region_grid=True,
                           save=True,
                           highlight1=[(5,124289158)],
                           pinpoint2=[(2,653874)]
                           )

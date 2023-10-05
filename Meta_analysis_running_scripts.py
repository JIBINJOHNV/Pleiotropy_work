
import os

####Munge sumstat
################################################---------------------------Summary sat synchronisation using Mungesumstat------------------------------------------------------#########################################
##PGC_SCZ_EUR
os.system('''zgrep -v "##" /home/jjohn41/Analysis/MetaAnalysis/GWAS_Raw_Data/PGC3_SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz  \
            > /home/jjohn41/Analysis/MetaAnalysis/GWAS_Raw_Data/PGC3_SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv''')

os.system('''nohup python3 ~/Softwares/CustomScripts/MungeSumstats_Running_Pysparq.py \
--inputfile /home/jjohn41/Analysis/MetaAnalysis/GWAS_Raw_Data/PGC3_SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv \
--outputfile SCZ_PGC3_Primary_22_EUR \
--infofilter 0.3 \
--frequencyfilter 0.001 \
--traitsnamee SCZ_PGC3_Primary_22_EUR \
--rsidcolumn ID \
--chr CHROM --position POS \
--a1 A1 --a2 A2 \
--allefreq FCON --info IMPINFO \
--beta BETA --se SE --pvalue PVAL \
--casescolumns  NCAS --controlscolumns NCON > PGC3_SCZ_wave3.european.autosome.public.v3.vcf.out 2>&1 &''')

os.system('''rm /home/jjohn41/Analysis/MetaAnalysis/GWAS_Raw_Data/PGC3_SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv''')

#Lam_et_al_2021_CognitiveTaskPerformance
os.system('''nohup python3 ~/Softwares/CustomScripts/MungeSumstats_Running_Pysparq.py \
--inputfile  /home/jjohn41/Analysis/MetaAnalysis/GWAS_Raw_Data/CTP_Lam_2021_EUR/Lam_et_al_2021_CognitiveTaskPerformance.tsv.gz \
--outputfile CTP_LAM_21_EUR \
--infofilter 0.00001 \
--frequencyfilter 0.00001 \
--traitsnamee CTP_LAM_21_EUR \
--rsidcolumn SNP \
--chr CHR --position BP \
--a1 A1 --a2 A2 \
--beta BETA --se SE \
--pvalue P --totalsamples 373617 > CTP_LAM_21_EUR_Lam_et_al_2021_CognitiveTaskPerformance.tsv.out 2>&1 &''')

#http://www.haplotype-reference-consortium.org/site
os.system('''bcftools annotate -a /home/jjohn41/Softwares/Resourses/dbsnps/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz \
                 -c INFO/AF CTP_LAM_21_EUR_MungeSumstats_Out_common.vcf.gz | bgzip -c > Annotated_output.vcf.gz''')

os.system('''tabix -p vcf Annotated_output.vcf.gz''')
os.system('''gatk VariantsToTable -F CHROM -F POS -F REF -F ALT -F AF -ASGF SNP -ASGF P -ASGF Z -ASGF N -V Annotated_output.vcf.gz -O MasterTable.tsv''')
os.system('''sed -i 's/CTP_LAM_21_EUR.//g' MasterTable.tsv''')

#cogent.iq.meta_1.tbl.snp.rsid.chr.bp.fuma.in.2.gz   ;COG_Lam_17_EUR
os.system('''nohup python3 ~/Softwares/CustomScripts/MungeSumstats_Running_Pysparq.py \
--inputfile MasterTable.tsv \
--outputfile CTP_LAM_21_EUR \
--infofilter 0.3 \
--frequencyfilter 0.001 \
--traitsnamee CTP_LAM_21_EUR \
--rsidcolumn SNP \
--chr CHROM --position POS \
--a1 ALT --a2 REF \
--zscore  Z --pvalue P \
--allefreq AF \
--samplescolumns N > CTP_LAM_21_EUR_Lam_et_al_2021_CognitiveTaskPerformance.tsv_2.out 2>&1 &''')


#GWAS_EA_excl23andMe.txt   ;EDU_Ex23ME_LEE_18_EUR
os.system('''nohup python3 ~/Softwares/CustomScripts/MungeSumstats_Running_Pysparq.py \
--inputfile /home/jjohn41/Analysis/MetaAnalysis/GWAS_Raw_Data/EDU_Ex23ME_LEE_18_EUR/GWAS_EA_excl23andMe.txt \
--outputfile EDU_Ex23ME_LEE_18_EUR \
--infofilter 0.3 \
--frequencyfilter 0.001 \
--traitsnamee EDU_Ex23ME_LEE_18_EUR \
--rsidcolumn MarkerName \
--chr CHR --position POS \
--a1 A1 --a2 A2 \
--beta Beta --se SE --pvalue Pval \
--allefreq EAF \
--totalsamples 766345 > EDU_Ex23ME_LEE_18_EUR_GWAS_EA_excl23andMe.txt.out 2>&1 &''')

os.system('''bcftools merge \
    SCZ_PGC3_Primary_22_EUR_MungeSumstats_Out_common.vcf.gz CTP_LAM_21_EUR_MungeSumstats_Out_common.vcf.gz \
    EDU_Ex23ME_LEE_18_EUR_MungeSumstats_Out_common.vcf.gz | bgzip -c > Recent_SCZ_COG_EDU.vcf.gz''')

os.system('''tabix -f -p vcf Recent_SCZ_COG_EDU.vcf.gz''')


#Pleio
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/1_PLEIO_filepreparation_pysparq.py \
        -sumstatvcffile /home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Input_File/Mungesumstat_Output/Recent_SCZ_COG_EDU.vcf.gz \
        -binarytraits SCZ_PGC3_Primary_22_EUR >  PLEIO_Recent_SCZ_COG_EDU_PLEIO_InputFilePreparation.out 2>&1 &''')

conda activate py27

os.system(f'''nohup python2 /home/jjohn41/Softwares/CustomScripts/2_PILEIO_Running.py -inputfile Pleio_input_list.txt >  PLEIO_Recent_SCZ_COG_EDU_LDSC.out 2>&1 &''')

conda activate base

os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/3_PILEIO_Running.py >  PLEIO_Recent_SCZ_COG_EDU_PLEIO_Analysis.out 2>&1 &''')


######--------------------------------Identify concordant, Discordant, SCZ spcific, EDU  Specific, COG specific---------------------------------------------------------######################
import os
import warnings
import argparse
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, when, array, arrays_zip,udf
from pyspark.sql.types import StringType, DateType, FloatType, DoubleType
from pyspark.sql.functions import concat_ws, array_join,array, arrays_zip, when, col,flatten,round,regexp_replace

effect="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/PLEIO_Results/PLEIO_META_output.blup.gz"
pleiovcf="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/PLEIO_Results/PLEIO_META_PLEIO_LSP_MungeSumstats_Out_common.vcf.gz"
pleiobetacutof=0.001
####vcf file
os.system(f"gatk VariantsToTable  -F CHROM -F POS -F REF -F ALT  -ASGF SNP -ASGF P  -V {pleiovcf} -O PLEIO_Effectsize_MasterTable.tsv")
spark = SparkSession.builder.config("spark.executor.memory", "30g").getOrCreate()

################--------------------------------------------------------------------------------Pleio out put---------------------------------------------------------------------------
vcfdf = spark.read.option("delimiter", "\t").option("header", True).csv("PLEIO_Effectsize_MasterTable.tsv")
new_cols=[x.replace(".","_") for x in vcfdf.columns ]
vcfdf = vcfdf.toDF(*new_cols)

vcfdf = vcfdf.drop("PLEIO_META_LSP_SNP", "PLEIO_META_LSP_P")
vcfdf = vcfdf.withColumn("PLEIO_META_PLEIO_P", when(vcfdf["PLEIO_META_PLEIO_P"] >=1, 0.999).otherwise(vcfdf["PLEIO_META_PLEIO_P"]))

effectdf = spark.read.option("delimiter", "\t").option("header", True).csv(effect)

# performing inner join on the two dataframes on the given columns
pleioout = vcfdf.join(effectdf, (col("PLEIO_META_PLEIO_SNP") == col("SNP")), "inner").drop("PLEIO_META_PLEIO_SNP")

Renamedict={"SCZ_PGC3_Primary_22_EUR":"SCZ_PLEIO_Beta","CTP_LAM_21_EUR":"COG_PLEIO_Beta","EDU_Ex23ME_LEE_18_EUR":"EDU_PLEIO_Beta"}
for old_name,new_name in Renamedict.items():
    pleioout = pleioout.withColumnRenamed(old_name,new_name)
    pleioout = pleioout.withColumn(new_name, col(new_name).cast("float"))
    pleioout = pleioout.withColumn(new_name, round(pleioout[new_name], 4)) 


##Converting the direction of SCHIZOPHRENI
pleioout=pleioout.withColumn("SCZ_PLEIO_Beta", col("SCZ_PLEIO_Beta")*-1)

#################-------------------------------------------------------Applying the directional effect and preparing input for FUMA ------------------------------------------------------------##############################
directionality_columns=['SCZ_PLEIO_Beta', 'COG_PLEIO_Beta', 'EDU_PLEIO_Beta']

positive_cols = arrays_zip(array(*[when(col(c) >=pleiobetacutof, c).otherwise("NA") for c in directionality_columns ]))
df_with_positive_cols = pleioout.withColumn("positive_cols", positive_cols)

negative_cols = arrays_zip(array(*[when(col(c) <=-pleiobetacutof, c).otherwise("NA") for c in directionality_columns]))
df_with_positive_cols = df_with_positive_cols.withColumn("negative_cols", negative_cols)


##Convert the  array to string
df_with_positive_cols = df_with_positive_cols \
    .withColumn("positive_cols_arr", flatten(array("positive_cols.0"))).withColumn("Pleio_positive", array_join("positive_cols_arr", ";"))

df_with_positive_cols = df_with_positive_cols \
    .withColumn("positive_cols_arr", flatten(array("negative_cols.0"))).withColumn("Pleio_negative", array_join("positive_cols_arr", ";"))

############Nonsignificant 
betacols1=[ x for x in df_with_positive_cols.columns if x.endswith("_Beta")]
betacols2=[ x for x in df_with_positive_cols.columns if x.endswith("_Beta")]

column_dict=dict(zip(betacols1,betacols2))

# Generate array of positive columns
positive_cols = arrays_zip(array(*[when((col(c) <pleiobetacutof) & (col(d) >-pleiobetacutof), c).otherwise("NA") for c, d in column_dict.items()]))
df_with_positive_cols = df_with_positive_cols.withColumn("positive_cols", positive_cols)

##Convert the  array to string
df_with_positive_cols = df_with_positive_cols \
    .withColumn("positive_cols_arr", flatten(array("positive_cols.0"))) \
    .withColumn("Pleio_NotSignificant", array_join("positive_cols_arr", ";"))

##Removing the unwantedcolumns
cols_to_drop = ["positive_cols_arr", "positive_cols","negative_cols"]
df_with_positive_cols=df_with_positive_cols.drop(*cols_to_drop)

##Removing the unwantedcolumns
df_with_positive_cols = df_with_positive_cols.select([regexp_replace(c, "NA;NA;NA", "NA").alias(c) for c in df_with_positive_cols.columns])
df_with_positive_cols = df_with_positive_cols.select([regexp_replace(c, "NA;NA", "NA").alias(c) for c in df_with_positive_cols.columns])

df_with_positive_cols = df_with_positive_cols.select([regexp_replace(c, ";NA", "").alias(c) for c in df_with_positive_cols.columns])
df_with_positive_cols = df_with_positive_cols.select([regexp_replace(c, "NA;", "").alias(c) for c in df_with_positive_cols.columns])

df_with_positive_cols = df_with_positive_cols.select([regexp_replace(c, "SCZ_PLEIO_Beta", "SCZ").alias(c) for c in df_with_positive_cols.columns])
df_with_positive_cols = df_with_positive_cols.select([regexp_replace(c, "COG_PLEIO_Beta", "COG").alias(c) for c in df_with_positive_cols.columns])
df_with_positive_cols = df_with_positive_cols.select([regexp_replace(c, "EDU_PLEIO_Beta", "EDU").alias(c) for c in df_with_positive_cols.columns])

df_with_positive_cols.coalesce(1).write.option("sep", "\t").csv('result.csv',header=True)
os.system(f"mv  result.csv/*csv AJHG_2019_Reanalysis_PLEIO_Output_with_positive_negative.tsv")
os.system("rm -rf result.csv")

spark.stop()

####Filtering convordant vs discordant
spark = SparkSession.builder.config("spark.executor.memory", "30g").getOrCreate()

################--------------------------------------------------------
df_with_positive_cols = spark.read.option("delimiter", "\t").option("header", True).csv("PLEIO_Output_with_positive_negative.tsv")


significant=df_with_positive_cols.filter(col("Pleio_NotSignificant")!="SCZ;COG;EDU")
significant_scz=significant.filter(col("Pleio_positive").contains("SCZ") | col("Pleio_negative").contains("SCZ"))

significant_scz_na=significant_scz.filter(col("Pleio_positive").contains("NA") | col("Pleio_negative").contains("NA"))

#significant_scz_na_noscz_na1=significant_scz_na.filter(col("Pleio_positive")!="SCZ" & col("Pleio_negative")!="SCZ")
significant_scz_na_noscz_na1=significant_scz_na.filter((col("Pleio_positive") != "SCZ") & (col("Pleio_negative") != "SCZ"))

significant_scz_na_noscz_na1_selected=significant_scz_na_noscz_na1.select(["CHROM","POS","REF","ALT","PLEIO_META_PLEIO_P","SNP","Pleio_positive","Pleio_negative","Pleio_NotSignificant"])

significant_scz_na_noscz_na1_selected.coalesce(1).write.option("sep", "\t").csv('result.csv',header=True)
os.system(f"mv  result.csv/*csv PLEIO_Output_with_Concordant_withPLeioBeta.tsv")
os.system("rm -rf result.csv")

###Discordant
dis_significant_scz=significant.filter(col("Pleio_positive").contains("SCZ") | col("Pleio_negative").contains("SCZ"))
dis_significant_scz2 = dis_significant_scz.filter((col("Pleio_positive") != "SCZ;COG;EDU") | (col("Pleio_negative") != "SCZ;COG;EDU"))
dis_significant_scz3 = dis_significant_scz2.filter((col("Pleio_positive") != "NA") & (col("Pleio_negative") != "NA"))
dis_significant_scz3_selected=dis_significant_scz3.select(["CHROM","POS","REF","ALT","PLEIO_META_PLEIO_P","SNP","Pleio_positive","Pleio_negative","Pleio_NotSignificant"])
dis_significant_scz3_selected.coalesce(1).write.option("sep", "\t").csv('result.csv',header=True)
os.system(f"mv  result.csv/*csv PLEIO_Output_with_Disconcordant_withPLeioBeta.tsv")
os.system("rm -rf result.csv")
dis_significant_scz3_selected.show()
significant_scz_na_noscz_na1_selected.groupby(['Pleio_positive', 'Pleio_negative', 'Pleio_NotSignificant']).count().show()
dis_significant_scz3_selected.groupby(['Pleio_positive', 'Pleio_negative', 'Pleio_NotSignificant']).count().show()


#############SCZ , EDU AND COG SPECIFIC
scz_specific=significant.filter( ((col("Pleio_positive")=="SCZ") & (col("Pleio_negative")=="NA")) | ((col("Pleio_negative")=="SCZ") & (col("Pleio_positive")=="NA")) )
scz_specific.coalesce(1).write.option("sep", "\t").csv('result.csv',header=True)
os.system(f"mv  result.csv/*csv PLEIO_Output_with_SCZ_Specific_withPLeioBeta.tsv")
os.system("rm -rf result.csv")

cog_specific=significant.filter( ((col("Pleio_positive")=="COG") & (col("Pleio_negative")=="NA")) | ((col("Pleio_negative")=="COG") & (col("Pleio_positive")=="NA")) )
cog_specific.coalesce(1).write.option("sep", "\t").csv('result.csv',header=True)
os.system(f"mv  result.csv/*csv PLEIO_Output_with_COG_Specific_withPLeioBeta.tsv")
os.system("rm -rf result.csv")

edu_specific=significant.filter( ((col("Pleio_positive")=="EDU") & (col("Pleio_negative")=="NA")) | ((col("Pleio_negative")=="EDU") & (col("Pleio_positive")=="NA")) )
edu_specific.coalesce(1).write.option("sep", "\t").csv('result.csv',header=True)
os.system(f"mv  result.csv/*csv PLEIO_Output_with_EDU_Specific_withPLeioBeta.tsv")
os.system("rm -rf result.csv")


cog_edu_together=significant.filter( ( (col("Pleio_positive").contains("COG") & (col("Pleio_positive").contains("EDU") )) |
                                       (col("Pleio_negative").contains("COG") & (col("Pleio_negative").contains("EDU") )) ))
cog_edu_together=cog_edu_together.filter( ~(col("Pleio_positive").contains("SCZ")) & ~(col("Pleio_negative").contains("SCZ")) )
cog_edu_together.coalesce(1).write.option("sep", "\t").csv('result.csv',header=True)
os.system(f"mv  result.csv/*csv PLEIO_Output_with_COG_EDU_together_Sameside_withPLeioBeta.tsv")
os.system("rm -rf result.csv")

cog_edu=significant.filter( ((col("Pleio_positive").contains("COG") | (col("Pleio_positive").contains("EDU") )) |
                    (col("Pleio_negative").contains("COG") | (col("Pleio_negative").contains("EDU") )) ))
cog_edu=cog_edu.filter( ~(col("Pleio_positive").contains("SCZ")) & ~(col("Pleio_negative").contains("SCZ")) )

cog_edu.coalesce(1).write.option("sep", "\t").csv('result.csv',header=True)
os.system(f"mv  result.csv/*csv PLEIO_Output_with_COG_or_EDU_withPLeioBeta.tsv")
os.system("rm -rf result.csv")

scz_specific.groupby(['Pleio_positive', 'Pleio_negative', 'Pleio_NotSignificant']).count().show()
cog_specific.groupby(['Pleio_positive', 'Pleio_negative', 'Pleio_NotSignificant']).count().show()
edu_specific.groupby(['Pleio_positive', 'Pleio_negative', 'Pleio_NotSignificant']).count().show()
cog_edu_together.groupby(['Pleio_positive', 'Pleio_negative', 'Pleio_NotSignificant']).count().show()
cog_edu.groupby(['Pleio_positive', 'Pleio_negative', 'Pleio_NotSignificant']).count().show()

significant_scz_na_noscz_na1_selected.groupby(['Pleio_positive', 'Pleio_negative', 'Pleio_NotSignificant']).count().show()
dis_significant_scz3.groupby(['Pleio_positive', 'Pleio_negative', 'Pleio_NotSignificant']).count().show()


#####################------------------------------------------------------FUMAINPUT FILE PRPARATION-----------------------------------------------------------------------##################

## Preparing the Directionality assigned files for FUMA
import pandas as pd
import glob,os

for file in glob.glob("*tsv"):
    df=pd.read_csv(file,sep="\t")
    df=df[['CHROM', 'POS', 'REF', 'ALT','SNP','PLEIO_META_PLEIO_P']].sort_values(by=["CHROM","POS"])
    columns=['CHROM', 'POS','REF',"EFFECT","SNP","P"]
    df.columns=['CHROM', 'POS','REF',"EFFECT","SNP","P"]
    df.to_csv(f"{file[:-4]}.txt",sep="\t",index=None)

##Preparing the PLEIO OUTPUT FOR FUMA analysis
import pandas as pd
import glob,os
os.system('gatk VariantsToTable  -F CHROM -F POS -F REF -F ALT  -ASGF SNP -ASGF P  -V PLEIO_META_PLEIO_LSP_MungeSumstats_Out_common.vcf.gz -O PLEIO_Effectsize_MasterTable.tsv')
df=df[['CHROM', 'POS', 'REF', 'ALT','PLEIO_META_PLEIO.SNP', 'PLEIO_META_PLEIO.P']].sort_values(by=["CHROM","POS"])
columns=['CHROM', 'POS','REF',"EFFECT","SNP","P"]
df.columns=['CHROM', 'POS','REF',"EFFECT","SNP","P"]
df.to_csv(f"COG_EDU_SCZ_PLEEIO_Analysis_Results.txt",sep="\t",index=None)

###Prepare FUMA input files from original SCZ, CTP  and edu GWAS
import pandas as pd
import glob,os
import numpy as np

os.system('gatk VariantsToTable -F CHROM -F POS -F REF -F ALT -ASGF FRQ -ASGF SNP -ASGF INFO -ASGF P -ASGF Z -V Recent_SCZ_COG_EDU.vcf.gz -O MasterTable.tsv')
df=pd.read_csv("MasterTable.tsv",sep="\t")
df['MAF']=np.where(df['SCZ_PGC3_Primary_22_EUR.FRQ']>0.05,1-df['SCZ_PGC3_Primary_22_EUR.FRQ'],df['SCZ_PGC3_Primary_22_EUR.FRQ'])
df=df[(df['SCZ_PGC3_Primary_22_EUR.INFO']>0.3) & (df['MAF']>0.001)]

columns=['CHROM', 'POS','REF',"EFFECT","SNP","P","Z"]

CTP=df[['CHROM', 'POS', 'REF', 'ALT','CTP_LAM_21_EUR.SNP', 'CTP_LAM_21_EUR.P', 'CTP_LAM_21_EUR.Z']]
CTP=CTP[~CTP['CTP_LAM_21_EUR.P'].isna()]
CTP.columns=columns
CTP.to_csv(f"CTP_GWAS_Results.txt",sep="\t",index=None)


EDU=df[['CHROM', 'POS', 'REF', 'ALT','EDU_Ex23ME_LEE_18_EUR.SNP','EDU_Ex23ME_LEE_18_EUR.P', 'EDU_Ex23ME_LEE_18_EUR.Z']]
EDU=EDU[~EDU['EDU_Ex23ME_LEE_18_EUR.P'].isna()]
EDU.columns=columns
EDU.to_csv(f"EDU_Results.txt",sep="\t",index=None)


SCZ=df[['CHROM', 'POS', 'REF', 'ALT','SCZ_PGC3_Primary_22_EUR.SNP','SCZ_PGC3_Primary_22_EUR.P', 'SCZ_PGC3_Primary_22_EUR.Z']]
SCZ=SCZ[~SCZ['SCZ_PGC3_Primary_22_EUR.P'].isna()]
SCZ.columns=columns
SCZ.to_csv(f"PGC3_SCZ_Results.txt",sep="\t",index=None)




#####################################-------------------------FUMA OUTPUT  LOCI MERGINGanndideentifying common LOCI -----------------------------###################################################
import glob
import pandas as pd
from pybedtools import BedTool
import numpy as np

files=glob.glob("/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/FUMA_Result/Concordant_Discordant/*GenomicRiskLoci.txt")

MainDf_list=list()
MainDf=pd.DataFrame()

for file in files:
    df1=pd.read_csv(file,sep="\t")
    df1a=df1[['chr','start', 'end']]
    MainDf_list.append(df1a)

MainDf=pd.concat(MainDf_list)
MainDf=MainDf.sort_values(by=["chr","start","end"]).drop_duplicates()
MainDfbed = BedTool.from_dataframe(MainDf)
MainDf_merged=MainDfbed.merge()
FinalResulktsdf=MainDf_merged.to_dataframe().copy()

for file in files:
    df1=pd.read_csv(file,sep="\t")
    colprefix=file.split("/")[-1].split("_")[2]+":"
    df1a=df1[['chr','start', 'end','nIndSigSNPs','nLeadSNPs','rsID','p']]
    df1abed = BedTool.from_dataframe(df1a)
    df1bed_merged=MainDf_merged.intersect(df1abed, wa=True, wb=True,loj=True)
    resultdf=df1bed_merged.to_dataframe()
    resultdf.columns=["chrom","start","end"]+[colprefix+x for x in df1a.columns  ]
    FinalResulktsdf=pd.merge(FinalResulktsdf,resultdf,on=['chrom', 'start', 'end'],how="outer")
    FinalResulktsdf=FinalResulktsdf.drop_duplicates()

df2=FinalResulktsdf.copy()
df2=df2.astype("str")

df3=df2.groupby(["chrom","start","end"]).agg({x:",".join for x in df2.iloc[:,3:].columns}).reset_index()

for Col in list(df3.iloc[:,3:].columns):
    df3[Col]=df3[Col].apply(lambda x: ";".join(set(x.split(","))))

df3=df3.sort_values(by=["chrom","start","end"]).drop_duplicates()
df3=df3.astype("str")
df3=df3.replace(".","NA")
df3=df3.replace("-1","NA")


condition=[ (df3['Concordant:chr']!="NA") & (df3['Discordant:chr']!="NA") ,
            ((df3['Concordant:chr']=="NA") & (df3['Discordant:chr']!="NA")),
            ((df3['Concordant:chr']!="NA") & (df3['Discordant:chr']=="NA")) ]

choicelist=["Both","Discordant Specific","Concordant Specific"]
df3['Categary']=np.select(condition,choicelist)
df3.to_csv("Fuma_Genomic_Loci_Comparison.csv",index=None)


Both=df3[df3['Categary']=="Both"]
Concordant=df3[df3['Categary']=="Concordant Specific"]
Discordant=df3[df3['Categary']=="Discordant Specific"]

Both.to_csv("Common_Fuma_Genomic_Loci.csv",index=None)
Concordant.to_csv("Concordant_Fuma_Genomic_Loci.csv",index=None)
Discordant.to_csv("Discordant_Fuma_Genomic_Loci.csv",index=None)




df3=pd.read_csv("Fuma_Genomic_Loci_Comparison.csv")
df3=df3.fillna("NA")
df5=df3[['chrom', 'start', 'end']+[x for x in df3.columns if x.endswith(":p")]]

for col in [x for x in df3.columns if x.endswith(":p")]:
    df5[col]=np.where(df5[col]=="NA",0,1)

df5[[ x for x in df5.columns if x.endswith(":p")]].to_csv("For_UpSetR_Analysis.csv",index=None)

df5['Common']=np.where( ( (df5['Concordant:p']==1)  & (df5['Discordant:p']==1)  ),1,0)
df5['Concordant:p']=np.where(df5['Common']==1,0,df5['Concordant:p'])
df5['Discordant:p']=np.where(df5['Common']==1,0,df5['Discordant:p'])

value_counts = {}
for column in df5.iloc[:,3:].columns:
    value_counts[column] = df5[column].value_counts()

grouped=pd.DataFrame(df5.groupby(['Common','Concordant:p','Discordant:p' ])['start'].count()).reset_index()

###Totalnumber ofloci 
loci_couts={}
for file in files:
    df1=pd.read_csv(file,sep="\t")
    loci_couts[file.split("/")[-1].split("_")[2]]=df1.shape[0]



####Vendiagram
df5[['chrom', 'start', 'end']]=df5[['chrom', 'start', 'end']].astype("str")
df5['value']=df5['chrom']+"_"+df5['start']+"_"+df5['end']
df5['Concordant:p']=np.where(df5['Concordant:p']==1,df5['value'],0)
df5['Discordant:p']=np.where(df5['Discordant:p']==1,df5['value'],0)
df5['Common']=np.where(df5['Common']==1,df5['value'],0)
df5=df5.replace(0,np.nan)
df5.to_csv("For_Vendiagram_Analysis.csv",index=None)



#####################################-------------------------MAGMA Analysis--------------------------------------------##############################

#######-----------------Prepare input for MAGMA Analysis--------------------------------------------###########
import pandas as pd

locifile="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/FUMA_Result/Concordant_Discordant/Fuma_Genomic_Loci_Comparison.csv"
concordant="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Concordant_Discordant_Assignment/PLEIO_Output_with_Concordant_withPLeioBeta.tsv"
discodant="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Concordant_Discordant_Assignment/PLEIO_Output_with_Disconcordant_withPLeioBeta.tsv"


locidf=pd.read_csv(locifile)
commdon_loci_df=locidf[locidf['Categary']=="Both"][['chrom', 'start', 'end']]

concordant_df=pd.read_csv(concordant,"\t")
concordant_df[['CHROM','POS']]=concordant_df[['CHROM','POS']].astype("int")
discodant_df=pd.read_csv(discodant,"\t")
discodant_df[['CHROM','POS']]=discodant_df[['CHROM','POS']].astype("int")


concordant_common=pd.DataFrame()

concordant_specific=concordant_df.copy()
discordant_specific=discodant_df.copy()


for index in range(commdon_loci_df.shape[0]):
    chr=commdon_loci_df.iloc[index,0]
    pos1=commdon_loci_df.iloc[index,1]
    pos2=commdon_loci_df.iloc[index,2]
    df2=concordant_df[((concordant_df['CHROM']==chr) & (concordant_df['POS']>=pos1) & (concordant_df['POS']<=pos2))]
    concordant_common=pd.concat([concordant_common,df2])
    concordant_specific=concordant_specific[~((concordant_specific['CHROM']==chr) & (concordant_specific['POS']>=pos1) & (concordant_specific['POS']<=pos2))]




for index in range(commdon_loci_df.shape[0]):
    chr=commdon_loci_df.iloc[index,0]
    pos1=commdon_loci_df.iloc[index,1]
    pos2=commdon_loci_df.iloc[index,2]
    df2=discodant_df[((discodant_df['CHROM']==chr) & (discodant_df['POS']>=pos1) & (discodant_df['POS']<=pos2))]
    concordant_common=pd.concat([concordant_common,df2])
    discordant_specific=discordant_specific[~((discordant_specific['CHROM']==chr) & (discordant_specific['POS']>=pos1) & (discordant_specific['POS']<=pos2))]




concordant_specific.to_csv("PLEIO_Output_with_Concordant_Specific_withPLeioBeta.tsv",sep="\t",index=None)
discordant_specific.to_csv("PLEIO_Output_with_Discordat_Specific_withPLeioBeta.tsv",sep="\t",index=None)
concordant_common.to_csv("PLEIO_Output_with_Both_ConcordantDiscordant_Specific_withPLeioBeta.tsv",sep="\t",index=None)

###################------------------------Exteded MHC Magma Analysis
import pandas as pd
import os,glob

files=glob.glob("/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Inputs/NoMHC_Extended_25Mb_35Mb/*tsv")

for file in files:
    df=pd.read_csv(file,sep="\t")
    df=df[~( (df["CHROM"]==6)  & (df["POS"]>=25000000) & (df["POS"]<=35000000))]
    df.to_csv(file,sep="\t",index=None)


InputPath_1="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Inputs/NoMHC_Extended_25Mb_35Mb/"
InputPath_2="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Inputs/NoMHC_Extended_25Mb_35Mb/"
OutPath="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Outputs/NoMHC_Extended_25Mb_35Mb/"
geneset="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/MAGAM_Genesets/NoMHC_Extended_25Mb_35Mb/"

##Concordant specific
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_1}PLEIO_Output_with_Concordant_Specific_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Concordant_Specific_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoMHC_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Concordant_Specific_Magma.log 2>&1 &''')

#Disconcordant specific
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_1}PLEIO_Output_with_Discordat_Specific_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Disconcordant_Specific_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoMHC_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Discordat_Specific_Magma.log 2>&1 &''')

#Disconcordant & Concordent common
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_1}PLEIO_Output_with_Both_ConcordantDiscordant_Specific_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Both_ConcordantDiscordant_Specific_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoMHC_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Both_ConcordantDiscordant_Magma.log 2>&1 &''')


##Concordant
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_2}PLEIO_Output_with_Concordant_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Concordant_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoMHC_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Concordant_Magma.log 2>&1 &''')

##Disconcordant
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_2}PLEIO_Output_with_Disconcordant_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Disconcordant_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoMHC_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Disconcordant_Magma.log 2>&1 &''')



############################--------------------------------------------------__##########################################################################--------------------------------
import os

InputPath_1="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Inputs/"
InputPath_2="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Inputs/"
OutPath="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Outputs/"
geneset="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/MAGAM_Genesets/"

##Concordant specific
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_1}PLEIO_Output_with_Concordant_Specific_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Concordant_Specific_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoMHC_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Concordant_Specific_Magma.log 2>&1 &''')

#Disconcordant specific
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_1}PLEIO_Output_with_Discordat_Specific_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Disconcordant_Specific_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoMHC_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Discordat_Specific_Magma.log 2>&1 &''')

#Disconcordant & Concordent common
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_1}PLEIO_Output_with_Both_ConcordantDiscordant_Specific_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Both_ConcordantDiscordant_Specific_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoMHC_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Both_ConcordantDiscordant_Magma.log 2>&1 &''')


##Concordant
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_2}PLEIO_Output_with_Concordant_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Concordant_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoMHC_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Concordant_Magma.log 2>&1 &''')

##Disconcordant
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_2}PLEIO_Output_with_Disconcordant_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Disconcordant_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoMHC_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Disconcordant_Magma.log 2>&1 &''')




############################----------------------------------With MHC__##########################################################################--------------------------------
import os

InputPath_1="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Inputs/"
InputPath_2="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Inputs/"
OutPath="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/Magma_Analysis/Magma_Outputs/No_ChrXY_WithMHC/"
geneset="/home/jjohn41/Analysis/MetaAnalysis/Analysis/Recent_SCZ_COG_EDU/MAGAM_Genesets/No_ChrXY_WithMHC/"

##Concordant specific
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_1}PLEIO_Output_with_Concordant_Specific_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Concordant_Specific_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Concordant_Specific_Magma.log 2>&1 &''')

#Disconcordant specific
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_1}PLEIO_Output_with_Discordat_Specific_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Disconcordant_Specific_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Discordat_Specific_Magma.log 2>&1 &''')

#Disconcordant & Concordent common
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_1}PLEIO_Output_with_Both_ConcordantDiscordant_Specific_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Both_ConcordantDiscordant_Specific_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Both_ConcordantDiscordant_Magma.log 2>&1 &''')


##Concordant
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_2}PLEIO_Output_with_Concordant_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Concordant_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Concordant_Magma.log 2>&1 &''')

##Disconcordant
os.system(f'''nohup python3 /home/jjohn41/Softwares/CustomScripts/magma_analysis_parallel.py \
    --inputfile {InputPath_2}PLEIO_Output_with_Disconcordant_withPLeioBeta.tsv \
    --outputfile {OutPath}PLEIO_Output_with_Disconcordant_Magma --chr CHROM --rsidcolumn SNP --position POS --pvalue  PLEIO_META_PLEIO_P --totalsamples 1000000  \
    --genesetfile {geneset}msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoChrXY_Minimum10_Selected.gmt  \
    --generefLoc {geneset}NCBI37.3.gene_Name_NoChrXY.loc >  {OutPath}PLEIO_Output_with_Disconcordant_Magma.log 2>&1 &''')







#############Merging gene wise association results


import numpy as np
import pandas as pd

df=pd.read_csv("MAGA_Gene_wISE_Pvalue_ConcordantDiscordantBoth.tsv",sep="\t")
genes=df[['GENE','Conncordant_P','Discordant_P','Both_P']]

condition = [
    ((genes['Conncordant_P'] < 2.631578947368421e-06) & (genes['Discordant_P'] < 2.631578947368421e-06) | (genes['Both_P'] < 2.631578947368421e-06)),
    ((genes['Conncordant_P'] < 2.631578947368421e-06) & (genes['Discordant_P'] > 2.631578947368421e-06) & (genes['Both_P'] > 2.631578947368421e-06)),
    ((genes['Conncordant_P'] > 2.631578947368421e-06) & (genes['Discordant_P'] < 2.631578947368421e-06) & (genes['Both_P'] > 2.631578947368421e-06))
]

choice = ["both", "concordat", "discordat"]

genes['Categary']=np.select(condition, choice)
genes['Categary']=genes['Categary'].str.replace("0","NA")
genes=genes.drop(["con", "disco",  "both"],axis=1)


genes.groupby("Categary")['GENE'].count().reset_index()

genes.to_csv("MAGA_Gene_wISE_Pvalue_ConcordantDiscordantBoth_CategaryAssigned.tsv",sep="\t",index=None)
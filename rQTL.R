#rQTL

snps_file_name=paste("snp2.txt"); #have to move the column names one space over, might need to add an extra space in the text file
covariates_file_name = paste("cvrresid.txt");
nk_resid_top_file_name= paste("nk_bfit_top_resid.txt")

#Load genotype data
snpnk = SlicedData$new();
#snpnk$fileDelimiter = "\t";      # the TAB character
snpnk$fileOmitCharacters = "NA"; # denote missing values;
snpnk$fileSkipRows = 1;          # one row of column labels
snpnk$fileSkipColumns = 1;       # one column of row labels
snpnk$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snpnk$LoadFile(snps_file_name);

# Covariates file name
# Set to character() for no covariates

#covariates_file_name =character();
cvrtnk = SlicedData$new();
cvrtnk$fileDelimiter = "\t"; # the TAB character
cvrtnk$fileOmitCharacters = "NA"; # denote missing values;
cvrtnk$fileSkipRows = 1; # one row of column labels
cvrtnk$fileSkipColumns = 1; # one column of row labels
if(length(covariates_file_name)>0) {
  cvrtnk$LoadFile(covariates_file_name);
}


nktop = SlicedData$new();
nktop$fileOmitCharacters = "NA"; # denote missing values;
nktop$fileSkipRows = 1;          # one row of column labels
nktop$fileSkipColumns = 1;       # one column of row labels
nktop$fileSliceSize = 2000;      # read file in slices of 2,000 rows
nktop$LoadFile(nk_resid_top_file_name);

meq_nktop_nkcvrt= Matrix_eQTL_engine(
  snps = snpnk, #all snps
  gene = nktop, #residuals for info
  cvrt = cvrtnk,# nkresiduals
  output_file_name = tempfile(),
  pvOutputThreshold = 1e-10, 
  useModel = modelLINEAR, 
  errorCovariance = numeric(), 
  verbose = TRUE,
  pvalue.hist = "qqplot");


View(meq_nktop_nkcvrt$all$eqtls)
nktop_hits<-meq_nktop_nkcvrt$all$eqtls
write.table(nktop_hits,'06_09_15_nktop_hits_nkcvrt.txt',sep="\t")
#plot(meq_nktop_nkcvrt, pch = 16, cex = 0.7)

#SNP SYMBOL
nktop_snps<-as.character(meq_nktop_nkcvrt$all$eqtls$snps)
View(nktop_snps)
nktop_snps<-as.numeric(nktop_snps)
View(nktop_snps)
class(nktop_snps)
nktop_snps_symbol<-data2_add@gtdata@snpnames[nktop_snps]
View(nktop_snps_symbol)


# SNP DATA

nksnp_data<-as.numeric(data2_add[,nktop_snps])
View(nksnp_data)
write.table(nksnp_data,'06_09_15_nksnp_nkcvrt_data.txt',sep='\t')


# NKGENE SYMBOL: Symbol or Top 10% Symbol 
nk_gene<-meq_nktop_nkcvrt$all$eqtls$gene
nk_gene<-as.character(nk_gene)
nk_gene<-as.numeric(nk_gene)
View(nk_gene)
nktop_gene_symbol<-nk_bfit_symbol[nk_gene]
View(nktop_gene_symbol)
write.table(nktop_gene_symbol,'06_09_15_nkgenesymbol_nkcvrt.txt',sep="\t")

#GENE EXPRESSION
nk_exp<-Exp[nk_gene,]
View(nk_exp)
write.table(nk_exp,'06_09_15_nk_nkcvrt_exp.txt',sep="\t")

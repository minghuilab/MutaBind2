library(forestFloor)
library(randomForest)

#f<-read.csv(file='sinputfiles/skempi_single_mutation_mutmwtb1000_4191_addpred_modelbase_addlen_ShortVersion.txt',header = TRUE, sep = '\t')
#rownames(f) <- paste(f$PDBid,f$Partner1,f$Partner2,f$MutChain,f$Mutation_PDB,f$label_dataset,sep = '_')

#label <- 'ddG_exp_ave~ddVdw+ddPb+Foldx_fold+provean+ACC_wt+ACCmp_wt+len_mp_contact_op_wt10'
#data <- f[, c('ddG_exp_ave',strsplit(strsplit(gsub("\n","",label), '~')[[1]][2], '[+]')[[1]])]
# 必须要建立模型，不能直接load
#set.seed(100)
#mutabinds.rf <- randomForest(as.formula(label), data = data ,keep.inbag=TRUE,ntree = 500)
#save.image("sinputfiles/mutabinds.RData")


load("sinputfiles/mutabinds.RData")
test <- read.csv(file='/home/cheny/mutabind2/sexample/sexample.input.cleaned.outdata',header = TRUE, sep = '\t') 
test_info <- read.csv(file='/home/cheny/mutabind2/sexample/sexample.sunddg',header = TRUE, sep = '\t')
test_features <- test[,9:(ncol(test)-1)]
test_features$ACC_wt <- as.integer(test_features$ACC_wt)
test_features$ACCmp_wt <- as.integer(test_features$ACCmp_wt)

ff <- forestFloor(mutabinds.rf,data[-1],test_features,bootstrapFC = TRUE, calc_np = TRUE)
fc <- round(ff$FCmatrix,4)
fc <- as.data.frame(fc)
colnames(fc) <- c('DDE_vdw','DDG_solv','DDG_fold','CS','SA_com_wt','SA_part_wt','N_cont_wt','bootstrapFC')
fc$dE_vdw_wt <- 0
fc <- fc[,c('DDE_vdw','DDG_solv','DDG_fold','CS','SA_com_wt','SA_part_wt','N_cont_wt','dE_vdw_wt','bootstrapFC')]
#fc$Sum_Contribution <- rowSums(fc)

result <- cbind(test_info,fc[(nrow(fc)-nrow(test)+1):nrow(fc),])
result$Model_Bias <- 0.4286  # 0.42862
#result$MutaBindS <- result$Sum_Contribution+result$Model_Bias
write.table(result,file = '/home/cheny/mutabind2/sexample/sexample.input.cleaned.outdata.contribution',sep = '\t',col.names = TRUE,row.names = FALSE)
#write.table(result,file = '2019051700102656321324162/2019051700102656321324162.input.cleaned.outdata.final.contribution',sep = '\t',col.names = TRUE,row.names = TRUE)

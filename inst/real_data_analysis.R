library(tidyverse)
library(corrplot)
library(randomForestSRC)
library(ACAT)

data = rio::import(".../pca10.csv")
propen_score = data$propen_score
path_score = rio::import(".../FHS_ssGSEA.csv")
FHS_gene = rio::import(".../FHS_data_gene_match.csv")

df = path_score
rownames(df) <- df[, 1]
df = df[,-1]
df_t <- as.data.frame(t(df), stringsAsFactors = FALSE)
df_t$shareid <- rownames(df_t)
rownames(df_t) <- NULL
df_t <- df_t[, c("shareid", setdiff(names(df_t), "shareid"))]

id = as.vector(as.matrix((path_score[1,2:5619])))

df_t[,1] = as.vector(id)

id_vector = c(id)
df_t$shareid = id_vector

merged_df <- merge(FHS_gene, df_t, by = "shareid")

df_new = merged_df

var_vector = c(1,3,38,39,40,41,42,49,62,53:54,68:77,17738:18061)

df_subset = df_new[,var_vector]
df_subset = na.omit(df_subset)

df_new[, "lipid"] <- NA
df_new$lipid = df_new$`(KEGG) Lipid and atherosclerosis`

name = c("ZBTB40","LAMC1","DEDD","PTPN7","MAL","SATB2","SIDT1","HTT","IL4","CD74","HIST1H3H","CFB",
         "AKAP12","LPAL2","PSMB1","ITGB8","MS4A2","SCCPDH","C1orf198","LMAN2L","DCUN1D4","PDGFRA","MFAP3L",
         "HIST1H4A","ABCA13","FLJ43692","SQLE","PPP2CB","CTTN","MMP8",
         "CD9","GPR84","FAM19A2","SPRED1","CETP","BPI","TGM2","NGFRAP1")

df_subset = df_new %>% select(var_vector,name,"lipid")
df_subset = na.omit(df_subset)
data_new = df_subset
data_new$Sex = ifelse(data_new$Sex == 2, yes = 1, no = 0)


data0 = data_new[data_new$Sex==0,]
data1 = data_new[data_new$Sex==1,]


fit_m_trt0 <- rfsrc(as.formula(formula_str_M), data = data0, ntree = 5000)
fit_m_trt1 <- rfsrc(as.formula(formula_str_M), data = data1, ntree = 5000)


formula_str_MY <- paste("HDL9~",paste(covariate_names_mediator, collapse = " + "),"+idtype+Age8+BMI8+Eversmk8+Currsmk8+Alcohol8_1+Alcohol8_2+lipid")  #paste("HDL9 ~Sex+ID_2991860 + ID_3137120 + ID_3332276+idtype+Age8+Eversmk8+Currsmk8+Alcohol8_1+Alcohol8_2+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")

rf_y_1 = rfsrc(as.formula(formula_str_MY),data = data1, ntree = 5000, var.used = "all.trees")

n = dim(data_new)[1]
n_mediator = length(name)#
df_pred_M_trt_1 = predict(object=fit_m_trt1, newdata = data_new)
pred_M_trt_1 = matrix(rnorm(n * n_mediator), n, n_mediator)
for (i in 1:n_mediator){
  pred_M_trt_1[,i] = df_pred_M_trt_1$regrOutput[i][[1]][1]$predicted
}

n = dim(data_new)[1]
df_pred_M_trt_0 = predict(object=fit_m_trt0, newdata = data_new)
pred_M_trt_0 = matrix(rnorm(n * n_mediator), n, n_mediator)

for (i in 1:n_mediator){
  pred_M_trt_0[,i] = df_pred_M_trt_0$regrOutput[i][[1]][1]$predicted

}


data_trt1_m1 = data_new
data_trt1_m1[,name] = pred_M_trt_1
Y_trt1_m1 = predict(object = rf_y_1, newdata = data_trt1_m1)

data_trt1_m0 = data_new
data_trt1_m0[,name] = pred_M_trt_0
Y_trt1_m0 = predict(object = rf_y_1, newdata = data_trt1_m0)

indir_trt_effect = Y_trt1_m1$predicted - Y_trt1_m0$predicted
#data
data_new$indir_trt_effect = indir_trt_effect

dis3 = as.data.frame(matrix(rep(0,n*n),nrow = n))
for (i in 1:n){
  dis3[i,] = rep(as.numeric(data_new[i,"indir_trt_effect"]),n) - as.numeric(data_new$indir_trt_effect)
}
dis3 = (dis3)**2
tsne_result = Rtsne::Rtsne(dis3,dims=2,is_distance=TRUE,verbose=FALSE,max_iter = 4000, theta = 0)


rio::export(tsne_result$Y, ".../FHS_data_clustering_tsne.csv")


kmeans_result_2 = kmeans(tsne_result$Y, centers = 2, iter.max = 100, nstart = 10)
kmeans_result_3 = kmeans(tsne_result$Y, centers = 3, iter.max = 100, nstart = 10)
kmeans_result_4 = kmeans(tsne_result$Y, centers = 4, iter.max = 100, nstart = 10)
kmeans_result_5 = kmeans(tsne_result$Y, centers = 5, iter.max = 100, nstart = 10)
kmeans_result_6 = kmeans(tsne_result$Y, centers = 6, iter.max = 100, nstart = 10)
kmeans_result_7 = kmeans(tsne_result$Y, centers = 7, iter.max = 100, nstart = 10)

data_new$kmeans2 = as.factor(kmeans_result_2$cluster)
data_new$kmeans3 = as.factor(kmeans_result_3$cluster)
data_new$kmeans4 = as.factor(kmeans_result_4$cluster)
data_new$kmeans5 = as.factor(kmeans_result_5$cluster)
data_new$kmeans6 = as.factor(kmeans_result_6$cluster)
data_new$kmeans7 = as.factor(kmeans_result_7$cluster)

pval_flag_M = 1
for (j in 2:7){
  x_name = c("idtype","Age8","BMI8","Eversmk8","Currsmk8","Alcohol8_1","Alcohol8_2","lipid")#colnames(data_new)[c(1,3:6,9:20)]#paste0("Cov",1:100)
  c_name = paste0("kmeans",j)  #"kmeans3"
  profiles <- tree_fit(Y=data_new[,c_name], X=data_new[,x_name], seed=1234,maxdepth = 3)
  #length(profiles$trees)
  plot_profile(profiles$trees[[1]])
  if(length(profiles$trees)>0){
    for(l in 1:length(length(profiles$trees))){
      data_tmp = predict_path(profiles$trees[[l]], newdata = data_new)
      fit_lm1 = lm(HDL9~leaf+Sex+leaf*Sex, data = data_tmp)
      fit_lm0 = lm(HDL9~leaf+Sex, data = data_tmp)
      pvalY <- as.numeric(na.omit(stats::anova(fit_lm0, fit_lm1)[[6]]))

      med_p_list = NULL
      for (k in 1:length(name)){
        formu_m1 = paste0((name[k]), "~leaf+Sex+leaf*Sex")
        formu_m0 = paste0((name[k]), "~leaf+Sex")
        fit_lm1 = lm(as.formula(formu_m1), data = data_tmp)
        fit_lm0 = lm(as.formula(formu_m0), data = data_tmp)
        pvalM <- as.numeric(na.omit(stats::anova(fit_lm0, fit_lm1)[[6]]))
        med_p_list = c(med_p_list, pvalM)
      }
      weight = rep(1,length(name))
      p_acat = ACAT(med_p_list, weights = weight)

      if(p_acat<pval_flag_M){
        pval_flag_M = p_acat
        pval_flag_Y = pvalY
        #data_flag = data_tmp
        n_cluster_M = j
        profile_index_M = 1
        profiles_save_M = profiles
        profiles_tree_index = l
      }

    }
  }

}


rio::export(data_new, ".../FHS_data_clustering.csv")

profiles_save_M = rio::import(".../lipid.rds")
rio::export(profiles_save_M,".../lipid.rds")
plot_profile(profiles_save_M$trees[[profiles_tree_index]])

data_tmp1 = predict_path(profiles_save_M$trees[[profiles_tree_index]], newdata = data_new)
x_name = c("idtype","Age8","BMI8","Eversmk8","Currsmk8","Alcohol8_1","Alcohol8_2","lipid")#colnames(data_new)[c(1,3:6,9:20)]#paste0("Cov",1:100)#colnames(data_new)[c(1,3:6,9:20)]#paste0("Cov",1:100)
c_name = "leaf"
profiles <- tree_fit(Y=data_tmp1[,c_name], X=data_tmp1[,x_name], seed=1234)
plot_profile(profiles$trees[[1]])


data_new = predict_path(profiles$trees[[1]], newdata = data_tmp)

data_subtype3 = df_new[df_new$BMI8>=26&df_new$lipid>=0.18,]


## Subtype 3
#data_subtype3 = data_subtype3[-which(is.na(out)),]
data_subtype3 = data_subtype3[-which(is.na(data_subtype3$HDL9)),]
out = data_subtype3$HDL9
med = data_subtype3[,name]
x = data_subtype3$Sex
res3 = Rsq.measure(p=1/2,Y=out,M=med,X=x,method = "ALL")

## Subtype 4
out = data_subtype4$HDL9
med = data_subtype4[,name]
x = data_subtype4$Sex
res4 = Rsq.measure(p=1/2,Y=out,M=med,X=x,method = "ALL")

## Subtype 6
out = data_subtype6$HDL9
med = data_subtype6[,name]
x = data_subtype6$Sex
res6 = Rsq.measure(p=1/2,Y=out,M=med,X=x,method = "ALL")

## Subtype 7
out = data_subtype7$HDL9
med = data_subtype7[,name]
x = data_subtype7$Sex
res7 = Rsq.measure(p=1/2,Y=out,M=med,X=x,method = "ALL")

# ## Subtype 9
# out = data_subtype9$HDL9
# med = data_subtype9[,name]
# x = data_subtype9$Sex
# res9 = Rsq.measure(p=1/2,Y=out,M=med,X=x,method = "ALL")



data_new1 = predict_path(profiles$trees[[1]], newdata = data_tmp1)
summary(data_new1$subtype)

x_name = colnames(data_new)[c(1,3:6,9:20)]#paste0("Cov",1:100)
c_name = "leaf"#paste0()  #"kmeans3"
profiles <- tree_fit(Y=data_new1[,c_name], X=data_new1[,x_name], seed=1234)
#length(profiles$trees)
plot_profile(profiles$trees[[1]])
data_new2 = predict_path(profiles$trees[[1]], newdata = data_new1)






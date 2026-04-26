##### Analyze single mediator Real data  ####
library(mediation)
df = rio::import(".../expe3_continuous/p.csv")
quantile(df$p1, 0.1)
#quantile(df$p2, 0.1)
quantile(df$p1, 0.05)
result = rio::import(".../single_mediator_result_new.rds")

res_list = list()
p_list = list()
for (i in 1:183){
  result_i = result[[i]]

  pval_flag_M = result_i[1]
  pval_flag_Y = result_i[2]
  n_cluster_M = result_i[3]
  profile_index_M = result_i[4]
  profiles_save_M = result_i[5]

  if(!is.null(profiles_save_M)){
    if(pval_flag_M<0.01521167){
      res_list = c(res_list, i)
      p_list = c(p_list,pval_flag_M)
    }
  }

}

##
##
j=27
i=res_list[[j]]
cat(i)
result_i = result[[i]]
profiles_save_M = result_i[5]
plot_profile(profiles_save_M[[1]]$trees[[1]])

data_new = rio::import(".../pca10.csv")

#library(mediation)
data_tmp = predict_path(profiles_save_M[[1]]$trees[[1]], newdata = data_new)

x_name = colnames(data_tmp)[c(1,3:6,9:20)] #paste0("Cov",1:100)
c_name = "leaf"#paste0("kmeans",j)  #"kmeans3"
profiles <- tree_fit(Y=data_tmp[,c_name], X=data_tmp[,x_name], seed=1234)
data_tmp = predict_path(profiles$trees[[1]], newdata = data_tmp)
plot_profile(profiles$trees[[1]])
summary(data_tmp$subtype)

data_sub = data_tmp[data_tmp$subtype%in%c(5),]
id = colnames(data_tmp)[20+i]
cat(id)
dat2[dat2$transcript_cluster_id == id,13]

formu_med = paste0(colnames(data_tmp)[20+i], "~Sex")
formu_out = paste0("HDL9~Sex+",colnames(data_tmp)[20+i])
med_model <- lm(as.formula(formu_med), data = data_sub)
out_model <- lm(as.formula(formu_out), data = data_sub)
med_out <- mediate(med_model, out_model, treat = "Sex", mediator = colnames(data_tmp)[20+i], sims = 1000)
summary(med_out)

## total
formu_med = paste0(colnames(data_tmp)[20+i], "~Sex")
formu_out = paste0("HDL9~Sex+",colnames(data_tmp)[20+i])
med_model <- lm(as.formula(formu_med), data = data_tmp)
out_model <- lm(as.formula(formu_out), data = data_tmp)
med_out <- mediate(med_model, out_model, treat = "Sex", mediator = colnames(data_tmp)[20+i], sims = 1000)
summary(med_out)


#res_list

med_name = NULL
for (i in 1:length(res_list)){
  j = res_list[[i]]
  med_name = c(med_name, colnames(data_tmp)[20+j])
}

rio::export(med_name,".../single_mediator_name.rds")


data_geneid = rio::import(".../FHS_Gene_TranscriptionID_annotation.RData")
dat1 = rio::import(".../FHS_HDL_Sex_Offsrping-Micro_Point_2024.RData")
dat2 = rio::import(".../FHS_Gene_TranscriptionID_annotation.RData")#FHS_HDL_Sex_Offsrping-RNAseq_Point_2024.RData

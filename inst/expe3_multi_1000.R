#### Expe 3 multiple mediators sample size =1000####
#####  Global #######
##### Experiment 3 sample size 1000 ######
##### Indirect experiment ########
library(randomForestSRC)
library(tidyverse)
library(rpart)
library(MASS)
library(rpart.plot)
library(ACAT)
library(foreach)
library(snowfall)

f_exp1 <- function(seed){

  set.seed(seed)

  n <- 1000
  p <- 50
  q <- 10
  X <- matrix(rnorm(n * p), n, p)

  W <- rbinom(n, 1, 0.5)
  Y = rep(0,n)

  # 1)
  #X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("Cov", 1:p)
  beta_mat = NULL
  # 2)
  if (is.null(beta_mat)) {
    beta_mat <- matrix(0, nrow = p, ncol = q)
    beta_mat[1, ] <- 0.5 + seq(-0.15,0.15,by=(0.3/9))     #
    beta_mat[2, ] <- 0.5 + seq(0.15,-0.15,by=-(0.3/9))   #
    beta_mat[3, ] <- 1  + seq(-0.15,0.15,by=(0.3/9)) #
    beta_mat[4, ] <- 1  + seq(0.15,-0.15,by=-(0.3/9))  #
  }

  # 3)
  mu <- X %*% beta_mat
  rho = 0.1
  q = 10
  sigma2 = 0.1
  # 4)
  Sigma <- sigma2 * ((1 - rho) * base::diag(q) + rho * matrix(1, q, q))

  # 5)
  E <- mvrnorm(n, mu = rep(0, q), Sigma = Sigma)
  #


  # 5)
  M <- mu + E
  colnames(M) <- paste0("M", 1:q)

  for (i in 1:n){
    if(W[i]==1){
      for(j in 1:q){
        M[i,j] = M[i,j]+W[i] + rnorm(1,0,0.01)
      }
    }else{
      for(j in 1:q){
        M[i,j] = M[i,j] + rnorm(1,0,0.01)
      }
    }
    Y[i] = 1 + 0.5*X[i,1]+0.5*X[i,2] +rnorm(1,0,0.01)+ 0.5*X[i,3] +0.5*X[i,4]+0.1*sum(M[i,])

  }
  data_X = as.data.frame(X)
  colnames(data_X) = paste0("Cov",1:50)
  #data$TRT = W
  data_M = M
  colnames(data_M) = paste0("M",1:q)
  data = cbind.data.frame(data_X,data_M)

  data$Y = Y
  data$propensity_score = p_treat
  data$TRT = W
  data$TRT = as.factor(data$TRT)

  covariate_names <- paste0("Cov", 1:50)
  mediator_names <-  paste0("M", 1:q)

  X = as.matrix(data[,1:50])
  y = as.factor(data$TRT)

  data0 = data[data$TRT==0,]
  data1 = data[data$TRT==1,]

  formula_str <- paste("Multivar(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10)", " ~TRT+", paste(covariate_names, collapse = " + "))

  fit_M_1 = rfsrc(as.formula(formula_str),
                  data = data1, ntree = 2000)
  fit_M_0 = rfsrc(as.formula(formula_str),
                  data = data0, ntree = 2000)
  formula_str_M <- paste("Y ~M1+M2+M3+M4+M5+M6+M7+M8+M9+M10+", paste(covariate_names, collapse = " + "))
  rf_y_1 = rfsrc(as.formula(formula_str_M),data = data1, ntree = 2000)

  n = dim(data)[1]
  n_mediator = q #
  df_pred_M_trt_1 = predict(object=fit_M_1, newdata = data)
  pred_M_trt_1 = matrix(rnorm(n * n_mediator), n, n_mediator)
  for (i in 1:n_mediator){
    pred_M_trt_1[,i] = df_pred_M_trt_1$regrOutput[i][[1]][1]$predicted
  }

  df_pred_M_trt_0 = predict(object=fit_M_0, newdata = data)
  pred_M_trt_0 = matrix(rnorm(n * n_mediator), n, n_mediator)

  for (i in 1:n_mediator){
    pred_M_trt_0[,i] = df_pred_M_trt_0$regrOutput[i][[1]][1]$predicted

  }

  data_trt1_m1 = data
  data_trt1_m1[,c(51:60)] = pred_M_trt_1#
  Y_trt1_m1 = predict(object = rf_y_1, newdata = data_trt1_m1)

  data_trt1_m0 = data#
  data_trt1_m0[,c(51:60)] = pred_M_trt_0#
  Y_trt1_m0 = predict(object = rf_y_1, newdata = data_trt1_m0)

  indir_trt_effect = Y_trt1_m1$predicted - Y_trt1_m0$predicted
  data$indir_trt_effect = indir_trt_effect

  dis3 = as.data.frame(matrix(rep(0,n*n),nrow = n))
  for (i in 1:n){
    dis3[i,] = rep(as.numeric(data[i,"indir_trt_effect"]),n) - as.numeric(data$indir_trt_effect)
  }
  dis3 = (dis3)**2
  tsne_result = Rtsne::Rtsne(dis3,dims=2,is_distance=TRUE,verbose=FALSE,max_iter = 4000, theta = 0)

  kmeans_result_2 = kmeans(tsne_result$Y, centers = 2, iter.max = 100, nstart = 10)
  kmeans_result_3 = kmeans(tsne_result$Y, centers = 3, iter.max = 100, nstart = 10)
  kmeans_result_4 = kmeans(tsne_result$Y, centers = 4, iter.max = 100, nstart = 10)
  kmeans_result_5 = kmeans(tsne_result$Y, centers = 5, iter.max = 100, nstart = 10)
  kmeans_result_6 = kmeans(tsne_result$Y, centers = 6, iter.max = 100, nstart = 10)
  kmeans_result_7 = kmeans(tsne_result$Y, centers = 7, iter.max = 100, nstart = 10)

  data$kmeans2 = as.factor(kmeans_result_2$cluster)
  data$kmeans3 = as.factor(kmeans_result_3$cluster)
  data$kmeans4 = as.factor(kmeans_result_4$cluster)
  data$kmeans5 = as.factor(kmeans_result_5$cluster)
  data$kmeans6 = as.factor(kmeans_result_6$cluster)
  data$kmeans7 = as.factor(kmeans_result_7$cluster)


  pval_flag_M = 1
  for (j in 2:7){
    x_name = paste0("Cov",1:50)#
    c_name = paste0("kmeans",j)  #
    profiles <- tree_fit(Y=data[,c_name], X=data[,x_name], seed=1234)

    if(length(profiles$trees)>0){
      for(l in 1:length(length(profiles$trees))){

        data_tmp = predict_path(profiles$trees[[l]], newdata = data)
        fit_lm1 = lm(Y~leaf+TRT+leaf*TRT, data = data_tmp)
        fit_lm0 = lm(Y~TRT+leaf, data = data_tmp)
        pvalY <- as.numeric(na.omit(stats::anova(fit_lm0, fit_lm1)[[6]]))

        med_p_list = NULL
        for (k in 1:q){

          formu_m1 = paste0("M",k, "~leaf+TRT+leaf*TRT")
          formu_m0 = paste0("M",k, "~TRT+leaf")
          fit_lm1 = lm(formu_m1, data = data_tmp)
          fit_lm0 = lm(formu_m0, data = data_tmp)
          pvalM <- as.numeric(na.omit(stats::anova(fit_lm0, fit_lm1)[[6]]))
          med_p_list = c(med_p_list, pvalM)
        }
        weight = rep(1,q)
        p_acat = ACAT(med_p_list, weights = weight)

        if(p_acat<pval_flag_M){
          pval_flag_M = p_acat
          pval_flag_Y = pvalY
          n_cluster_M = j
          profile_index_M = l
          profiles_save_M = profiles
        }

      }
    }

  }

  data_save_path = paste0(".../data_seed_",seed,".csv")
  rio::export(data, data_save_path)
  result = list()
  save_path = paste0(".../seed_",seed,".rds")
  result[[1]] = pval_flag_M
  result[[2]] = n_cluster_M
  result[[3]] = profile_index_M  #data_flag
  result[[4]] = profiles_save_M
  result[[5]] = pval_flag_Y

  names(result) = c("pval_M","n_cluster_M","profile_index_M","profiles_save_M","pval_flag_Y")
  rio::export(result, save_path)

}

library(Rtsne)
library(rio)
library(snowfall)
library(parallel)
library(MASS)
my_new_folder = ".../expe3_continuous_1000_multi"
if (!dir.exists(my_new_folder)) {
  dir.create(my_new_folder)
}

sfInit(parallel = TRUE, cpus = (detectCores() - 4))
sfLibrary(dplyr)
sfLibrary(rio)
sfLibrary(rpart)
sfLibrary(randomForestSRC)
sfLibrary(base)
sfLibrary(glmnet)
sfLibrary(MASS)
sfLibrary(Rtsne)
sfLibrary(ACAT)
sfLibrary(mhighlearner)

result_vector = sfLapply(1:100, f_exp1)


sfStop()





####### calculate coefficients in mediation models (FHS) #####

fhs_3=rio::import(".../fhs3.csv")
fhs_4=rio::import(".../fhs4.csv")
fhs_6=rio::import(".../fhs6.csv")
fhs_7=rio::import(".../fhs7.csv")
fhs_all = rbind.data.frame(fhs_3,fhs_4,fhs_6,fhs_7)

data_resample <- fhs_7 #

name = c("ZBTB40","LAMC1","DEDD","PTPN7","MAL","SATB2","SIDT1","HTT","IL4","CD74","HIST1H3H","CFB",
         "AKAP12","LPAL2","PSMB1","ITGB8","MS4A2","SCCPDH","C1orf198","LMAN2L","DCUN1D4","PDGFRA","MFAP3L",
         "HIST1H4A","ABCA13","FLJ43692","SQLE","PPP2CB","CTTN","MMP8",
         "CD9","GPR84","FAM19A2","SPRED1","CETP","BPI","TGM2","NGFRAP1")
M = data_resample[,name]

X = data_resample$Sex
Y = data_resample$HDL9



Y = scale(Y, center=T, scale=T)
X = scale(X, center=T, scale=T)
M = scale(M, center=T, scale=T)

tmp = data.frame(cbind(Y,X,M))
names(tmp)[1] = "Y"
names(tmp)[2] = "X"
fit.xm = lm(Y~X+M)
summary(fit.xm)

coef_tbl <- tidy(fit.xm)
write_xlsx(coef_tbl, ".../fhs_model_results_7.xlsx")


fit.mx = lm(M~X)
summary(fit.mx)
coef_tbl <- tidy(fit.mx)
coef_tbl <- coef_tbl[seq(2, nrow(coef_tbl), by = 2), ]
write_xlsx(coef_tbl, ".../fhs_model_results_7_mx.xlsx")
View(coef_tbl)

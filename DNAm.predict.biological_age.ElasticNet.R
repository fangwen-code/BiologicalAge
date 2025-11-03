library(limma)
library(glmnet)
library(caret)
library(tibble)
library(dplyr)

m_beta <- read.table("CpG.value.txt", sep="\t", header = TRUE, row.names = 1, check.names = FALSE)
m_pheno <- read.table("pheno.info.txt", sep="\t", header = TRUE, row.names = 1, check.names = FALSE)

common_samples <- intersect(colnames(m_beta), rownames(m_pheno))

m_beta <- m_beta[, common_samples, drop = FALSE]
m_pheno <- m_pheno[common_samples, , drop = FALSE]

### 筛选与衰老相关的位点
# ~0 表示不添加截距项，直接建模各变量的效应, 可以根据数据情况添加其他协变量，如细胞比例
# 细胞比例可通过ENmix包的estimateCellProp()函数进行计算
design <- model.matrix(~ 0 + Age + Gender2, data = m_pheno[, c("Age", "Gender2")])

# limma 线性回归与经验贝叶斯检验
fit <- lmFit(m_beta, design)
# 经验贝叶斯调整（提高小样本下的统计效能）
fit <- eBayes(fit)

# 提取年龄相关统计量并筛选DMPs
# 提取年龄（Age）的效应：包括系数（ΔBeta/year）、p值等
# logFC列即 ΔBeta/year（每增加 1 岁，beta 值的变化量），正值表示随年龄增加甲基化升高，负值表示降低。
age_result <- topTable(fit, coef = "Age", number = Inf, adjust.method = "BH")
age_result <- rownames_to_column(age_result, var = "CpG")

# 筛选符合条件的结果：
# BH-adjusted P-value < 0.01 and |ΔBeta/year| > 0.002
age_related_DMPs <- subset(age_result, adj.P.Val < 0.01 & abs(logFC) > 0.002)

# 保留结果到文件
write.table(age_related_DMPs, "age_related_DMPs.txt",sep = "\t", row.names = FALSE, quote = FALSE)

### 数据准备
# 提取 DMPs 的 CpG 位点名
dmp_cpgs <- age_related_DMPs$CpG
dmp_beta <- m_beta[dmp_cpgs, ,drop=FALSE]

com_samples <- intersect(colnames(dmp_beta), rownames(m_pheno))

age <- m_pheno[com_samples,]$Age
covariates <- m_pheno[com_samples,]$Gender2

# glmnet 要求 X 行=样本，列=特征
dmp_beta_t <- t(dmp_beta)

X <- cbind(dmp_beta_t, covariates)
y <- age

### 划分训练集与测试集(8:2)
set.seed(42)
# 保留年龄分布一致性，用 createDataPartition 更稳定, 训练集占比80%
train_index <- createDataPartition(y = y, p = 0.8, list = FALSE, times = 1)

# 生成训练集和测试集
X_train <- X[train_index, , drop = FALSE]
y_train <- y[train_index]                
X_test <- X[-train_index, , drop = FALSE]
y_test <- y[-train_index] 

# 划分情况统计
cat("训练集样本数：", nrow(X_train), "  测试集样本数：", nrow(X_test), "\n")
cat("训练集年龄均值：", round(mean(y_train), 2), "  测试集年龄均值：", round(mean(y_test), 2), "\n")

# 标准化：用训练集的统计量（均值、标准差）标准化测试集
train_mean <- colMeans(X_train, na.rm = TRUE)
train_sd <- apply(X_train, 2, sd, na.rm = TRUE)

# 处理标准差为0的特征（避免除以0，直接替换为1）
train_sd[train_sd == 0] <- 1

X_train_std <- (X_train - train_mean) / train_sd
X_test_std <- (X_test - train_mean) / train_sd

# 针对缺失值NA用训练集均值插补
X_train_std_imputed <- apply(X_train_std, 2, function(col) {
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

### 构建Elastic Net模型
## 筛选最优参数
alpha_grid <- seq(0.1, 0.9, by = 0.1) # 从0.1到0.9搜索，步长为0.1, 一般可以给0.5
nfold_use <- 10
cv_results <- list()

# 遍历每个alpha,进行10折交叉验证
# family 参数的选择取决于因变量（y）的类型,type.measure 的选择严格依赖于 family
# 连续数值选gaussian(默认值)，对应线性回归；type.measure可选mse, mae
# 二分类变量选binomial，对应逻辑回归；type.measure可选deviance, class, auc
# 多分类变量选multinomial，对应多分类逻辑回归；type.measure可选deviance, class
# 计数变量选poisson，type.measure可选deviance
for (alpha in alpha_grid) {
  cv_fit <- cv.glmnet(
    x = X_train_std_imputed,
    y = y_train,
    alpha = alpha,
    nfolds = nfold_use,
    family = "gaussian",
    type.measure = "mae",
    standardize = FALSE
  )
  cv_results[[as.character(alpha)]] <- cv_fit
}

# $cvm 存储每个 lambda 对应的评估指标，提取每个 alpha 对应的最小交叉验证 MAE
min_mae <- sapply(cv_results, function(x) min(x$cvm))

# 选择最小 MAE 对应的最优 alpha
best_alpha <- as.numeric(names(min_mae)[which.min(min_mae)])
best_cv_fit <- cv_results[[as.character(best_alpha)]]

# 选择最小 MAE 对应的最优 lambda
best_lambda <- best_cv_fit$lambda.min

# 输出最优超参数
cat("最优 alpha：", best_alpha, "\n")
cat("最优 lambda：", round(best_lambda, 6), "\n")
cat("训练集交叉验证最小 MSE：", round(min(best_cv_fit$cvm), 4), "\n")

# 基于最优参数构建最终模型
final_enet <- glmnet(
  x = X_train_std_imputed, 
  y = y_train, 
  alpha = best_alpha, 
  lambda = best_lambda, 
  family = "gaussian", 
  standardize = FALSE
)

## 模型评估
# 针对测试集缺失值NA用训练集的均值插补, 因X_test_std已经用训练集的均值和sd进行标准化了，所以，现在的缺失值用0进行插补
X_test_std_imputed <- X_test_std
for (j in 1:ncol(X_test_std_imputed)) {
  X_test_std_imputed[is.na(X_test_std_imputed[, j]), j] <- 0  # 标准化后均值为0，直接填0
}

y_pred_test <- predict(final_enet, newx = X_test_std_imputed)
y_pred_test <- as.vector(y_pred_test)                 

# R²：决定系数；MSE：均方误差；MAE：平均绝对误差
test_metrics <- data.frame(
  Metric = c("R²", "MSE", "MAE"),  
  Value = c(
    cor(y_test, y_pred_test)^2,             
    mean((y_test - y_pred_test)^2),           
    mean(abs(y_test - y_pred_test))           
  )
)

cat("测试集模型评估结果：\n")
print(test_metrics)

## 可视化
plot(
  x = y_test, 
  y = y_pred_test, 
  xlab = "Chronological Age (years)", 
  ylab = "DNAm Age (years)", 
  pch = 16, col = "black", cex = 0.5
)

abline(lm(y_pred_test ~ y_test), col = "steelblue", lwd = 2)

text(x = max(y_test)*0.32, y = max(y_pred_test)*0.93, 
     paste0("R² = ", round(test_metrics$Value[1], 3)), 
     col = "black", cex = 0.6)  # 添加 R² 标注
text(x = max(y_test)*0.32, y = max(y_pred_test)*0.95, 
     paste0("MAE = ", round(test_metrics$Value[3], 3)), 
     col = "black", cex = 0.6)  # 添加 MAE 标注

## 提取模型关键特征
# Elastic Net 通过 L1 正则化将无关特征的系数压缩至 0，非零系数对应的特征是模型的关键位点 / 协变量
# 提取最终模型的系数
enet_coef <- coef(final_enet)
coef_df <- data.frame(
  Feature = rownames(enet_coef)[-1],  # 排除截距项（第一行）
  Coefficient = as.vector(enet_coef)[-1],  # 排除截距项的系数
  stringsAsFactors = FALSE
)

# 按系数绝对值排序（系数绝对值越大，对预测的贡献越大）
key_features <- subset(coef_df, Coefficient != 0)
key_features <- key_features[order(abs(key_features$Coefficient), decreasing = TRUE), ]

# 输出关键特征数量和前10个特征
cat("模型中关键特征数量（非零系数）：", nrow(key_features), "\n")
cat("前10个关键特征（按系数绝对值排序）：\n")
print(head(key_features, 10))

write.table(key_features, "ENet_key_features.txt", sep = "\t", row.names = FALSE, quote = FALSE)
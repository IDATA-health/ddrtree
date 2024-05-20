library(survival)
library(cmprsk)
# 读取数据
data1 <- read.csv("C:\\Users\\WF\\Desktop\\DR\\ddrtree.csv")

# 将时间字符串转换为年月格式
data1$DNtime <- as.Date(data1$DNtime, format = "%Y-%m-%d")
data1$DRtime <- as.Date(data1$DRtime, format = "%Y-%m-%d")

data1$time <- as.numeric(difftime(data1$DRtime, data1$DNtime, units = "days"))*12 / 365

data <- data1[data1$time >= 0 , ]

fstatus <- ifelse(data$DR == 1, 1, 0)

covariates <- data.frame(sbp = data$sbp,
                         dbp = data$dbp,
                         HbA1c = data$HbA1c,
                         total_cholesterol = data$total_cholesterol,
                         triglycerides = data$triglycerides,
                         Creatinine = data$Creatinine,
                         ALT = data$ALT,
                         HDLC = data$HDLC)

# 进行竞争风险回归分析
fit <- crr(ftime = data$time,
           fstatus = fstatus,
           cov1 = covariates,failcode = 1)

# 计算多因素的累积发生概率
predictions <- matrix(NA, nrow = nrow(data1), ncol = 1)  # 创建存储预测结果的数组

for (i in 1:nrow(data1)) {
  patient_data <- data1[i, ]
  covariates <- data.frame(sbp = patient_data$sbp,
                           dbp = patient_data$dbp,
                           HbA1c = patient_data$HbA1c,
                           total_cholesterol = patient_data$total_cholesterol,
                           triglycerides = patient_data$triglycerides,
                           Creatinine = patient_data$Creatinine,
                           ALT = patient_data$ALT,
                           HDLC = patient_data$HDLC)
  
  # 使用 predict.crr 函数计算累积发生概率
  prediction <- predict(fit, cov1 = covariates)
  
  # 存储预测结果的最后一个时间点的发生概率
  predictions[i] <- prediction[length(prediction)]
}

# 输出预测结果数组
predictions



# 小数点后四位数格式化函数
format_decimal <- function(x) {
  format(x, digits = 4)
}

# 格式化后的预测结果
predictions_formatted <- apply(predictions, 2, format_decimal)

# 将预测结果和data1合并
data1_with_predictions <- cbind(data1, probability = predictions_formatted)

# 将数据保存为新的CSV文件
write.csv(data1_with_predictions, file = "new_data1.csv", row.names = FALSE)












# 加载包
library(readxl)
library(dplyr)
library(writexl)
library(tidyr)
library(openxlsx)

# 设置工作目录
setwd("D:/000weng_lab/data_R/growth_curve/20250420")

# 定义文件名
file_name <- "label.xlsx"
sheet_name <- "auto_label"

# 读取数据并处理
data_map <- read_excel("map_property.xlsx", sheet = "map")
data_property <- read_excel("map_property.xlsx", sheet = "property")
label_map <- pivot_longer(data_map, cols = -c("row"), names_to = "col", values_to = "plasmid")
label_map$well <- gsub(" ", "", paste(label_map$row, label_map$col))
label_map <- label_map %>% filter(!is.na(plasmid))
label_map$plasmid <- gsub(" ", "", paste("pTY", label_map$plasmid))
label <- merge(label_map, data_property, by = "plasmid")
label <- label %>% select(-row, -col)
label <- label %>%
    group_by(plasmid) %>%
    mutate(replicate = as.character(row_number())) %>%
    ungroup()

# 检查文件是否存在并加载/创建工作簿
if (file.exists(file_name)) {
    wb <- loadWorkbook(file_name)
    removeWorksheet(wb, sheet_name) # 删除旧工作表
} else {
    wb <- createWorkbook()
}
addWorksheet(wb, sheet_name)
writeData(wb, sheet = sheet_name, x = label) # 正确写入数据

# 保存工作簿（修正关键步骤）
saveWorkbook(wb, file_name, overwrite = TRUE)

print(paste("Label file_name has been created and saved as:", file_name))

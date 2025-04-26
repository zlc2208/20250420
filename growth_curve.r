library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(writexl)
library(randomcoloR)
library(gridExtra)
library(ggrepel)
library(colorspace)
library(drc)
library(ggh4x)

theme_set(theme_bw())
setwd("D:/000weng_lab/data_R/growth_curve/20250420") # setwd() 函数用于设置工作目录

if (file.exists("output/od_long.xlsx")) {
  print("od_long.xlsx exists, loading it...")
  od_long <- read_excel("output/od_long.xlsx") # 读取数据
} else {
  dir.create(path = "output") # 创建输出目录
  ## 给OD值加上标签并转化为长格式:od_long
  data_od <- read_excel("data.xlsx") # 读取数据
  label <- read_excel("label.xlsx", sheet = "label") # 读取标签
  # OD值减去blank均值
  OD_0 <- data_od %>%
    filter(data_od$well == "OD_0") %>%
    summarise(across(where(is.numeric), mean))
  real_od <- do.call("rbind", lapply(seq_len(nrow(data_od)), function(i) {
    data_od[i, -1] - OD_0
  }))
  real_od$well <- data_od$well # real_od没有第一列故添加第一列
  real_od <- real_od[!(real_od$well == "OD_0"), ]
  # 按 well 列合并两个数据框，即加上标签信息，再转化为长格式
  od_melt <- merge(real_od, label, by = "well")
  od_long <- pivot_longer(od_melt,
    cols = -c("well", "plasmid", "RNA7SK", "HEXIM1", "DHFR3", "DHFR12", "meaning", "replicate"),
    names_to = "time",
    values_to = "value",
    names_transform = list(time = as.numeric), # 转化为numeric格式
    values_transform = list(value = as.numeric)
  )
  od_long$time <- sapply(od_long$time, function(x) {
    x / 6
  })
  od_long$value[od_long$value < 0] <- 0
  write_xlsx(od_long, path = "output/od_long.xlsx") # 将数据写入Excel文件
}
#----绘制生长曲线图----
# 生成数据
if (TRUE) {
  od_long$plasmidid <- paste(od_long$plasmid, od_long$meaning)
  data_control <- od_long %>% filter(meaning %in% c("positive", "negative"))
  data_control <- subset(data_control, select = -c(RNA7SK, HEXIM1, DHFR3, DHFR12))
  data_figure <- od_long %>% filter(!meaning %in% c("positive", "negative"))
  # 绘图
  # 步骤1：生成主图的自动颜色
  main_colors <- c(
    `pTY120-1 double negative` = "#D1F7EF",
    `pTY121-5 double negative` = "#E2A8E5",
    `pTY121-1 double negative` = "#B36EEB",
    `pTY102-1 double negative` = "#F6E7FD",
    `pTY102-2 double negative` = "#73EDC8",
    `pTY102-4 double negative` = "#D5AB74",
    `pTY102-3 double negative` = "#EFE0BB",
    `pTY104-1 double negative` = "#F3D76E",
    `pTY126-1 double positive` = "#AFFBDD",
    `pTY126-3 double positive` = "#D6ED58",
    `pTY126-2 double positive` = "#89D0EC",
    `pTY126-4 double positive` = "#41CBD1",
    `pTY127-1 double positive` = "#625DBA",
    `pTY127-3 double positive` = "#A9C983",
    `pTY127-2 double positive` = "#92A3EF",
    `pTY127-4 double positive` = "#A046F7",
    `pTY128-1 double positive` = "#9FBAC8",
    `pTY128-3 double positive` = "#AABAE5",
    `pTY128-2 double positive` = "#E75CD4",
    `pTY128-4 double positive` = "#68D484",
    `pTY129-1 double positive` = "#EF8E6F",
    `pTY129-3 double positive` = "#91F162",
    `pTY129-2 double positive` = "#F1699F",
    `pTY129-4 double positive` = "#FCB4B9"
  )
  if (FALSE) {
    main_colors <- distinctColorPalette(length(unique(data_figure$plasmidid)), runTsne = TRUE)
    names(main_colors) <- unique(data_figure$plasmidid)
    dput(main_colors)
  }
  # 步骤2：手动设置对照组颜色
  control_colors <- c(
    "pGJJ588 positive" = "red",
    "pGJJ336 positive" = "orange",
    "pGJJ584 negative" = "blue"
  )
  # 步骤3：合并颜色，确保对照组颜色覆盖同名类别
  all_colors <- main_colors
  all_colors[names(control_colors)] <- control_colors
  p_1 <- tryCatch(
    {
      ggplot(
        data = data_figure,
        aes(x = time, y = value, group = well, color = plasmidid, alpha = replicate)
      ) +
        facet_grid(DHFR3 + DHFR12 ~ `RNA7SK` + `HEXIM1`) +
        geom_line() +
        scale_color_manual(values = all_colors) +
        theme(theme_classic())
    },
    error = function(e) {
      message(paste("!!!Error in p_1:", e$message))
      return(NULL)
    }
  )
  ggsave("output/01.png", plot = p_1)
  p_2 <- tryCatch(
    {
      p_1 + # 绘制分面对应的对照组曲线
        geom_line(
          data = data_control,
          aes(x = time, y = value, group = well, color = plasmidid, alpha = replicate),
          inherit.aes = FALSE,
          linewidth = 0.1, linetype = "dashed"
        )
    },
    error = function(e) {
      message(paste("!!!Error in p_2:", e$message))
      return(NULL)
    }
  )
  ggsave("output/02.png", plot = p_2)
}
#----线性拟合求生长速率----
if (TRUE) {
  if (file.exists("output/growth_rate.xlsx")) {
    print("growth_rate.xlsx exists, loading it...")
    gr_all <- read_excel("output/growth_rate.xlsx")
  } else {
    dir.create(path = "output", showWarnings = FALSE) # 创建输出目录
    gr_all <- do.call("rbind", lapply(unique(od_long$well), function(x) {
      ds <- od_long[od_long$well == x & od_long$time <= 75, ]
      max_od <- max(ds$value)
      mid_od <- max_od - ((max_od - min(ds$value)) / 2)
      time_odmax <- ds$time[ds$value == max_od]
      time_odmid <- max(ds$time[ds$value <= mid_od & ds$time <= time_odmax[1]]) # mid_od不一定存在，且存在衰退期使odmax不止存在一个值，故取第一个
      lmds <- data.frame(t = ds$time[ds$time <= time_odmid + 4 & ds$time >= time_odmid - 4]) # 从 ds 数据框中筛选出时间在 time_odmid 前后4h内的数据，并将这些时间点存储到新的数据框 lmds 中，作为 t 列。
      lmds$lgods <- log(ds$value[ds$time %in% lmds$t])
      lmfit <- lm(formula = lgods ~ t, data = lmds)
      lmds$fitted_data <- fitted(lmfit)
      r_squared <- summary(lmfit)$r.squared
      if (TRUE) { # 是否绘制详细拟合图
        p1 <- ggplot(lmds, aes(x = t, y = lgods)) +
          geom_point(color = "blue") + # 绘制原始数据点
          geom_line(aes(y = fitted_data), color = "red") + # 绘制拟合曲线
          annotate("text",
            x = max(lmds$t) - 2.5, y = max(lmds$fitted) - 0.04,
            label = paste(
              "y=", round(lmfit$coefficients[1], 4), "+", round(lmfit$coefficients[2], 4), "*x", "\n",
              "R² =", round(r_squared, 3), "\n",
              "time_odmid=", round(time_odmid, 1), "h", "\n",
              "growth_rate=", round(lmfit$coefficients[2], 4)
            ), hjust = 1
          ) + # 添加生长速率和R^2值
          labs(
            title = paste("Log(OD_600 Value)-Time_", x, ds$plasmid[1]),
            x = "Time/h",
            y = "Log(OD_600 Value)"
          )
        p2 <- ggplot(lmds, aes(x = t, y = ds$value[ds$time %in% lmds$t])) +
          geom_point(color = "blue") + # 绘制原始数据点
          labs(
            title = paste("OD_600 Value - Time_", x, ds$plasmid[1]),
            x = "Time/h",
            y = "OD_600 Value"
          )
        print(paste(x, ":"))
        combined_plot <- grid.arrange(p1, p2, ncol = 2)
        ggsave(paste("output/", x, ds$plasmidid[1], ".png"), combined_plot)
      } else {
        print("No plot")
      }
      data.frame(
        well = x, plasmid = ds$plasmid[1], replicate = ds$replicate[1], growth_rate = lmfit$coefficients[2], time_odmid = time_odmid,
        `RNA7SK` = ds$`RNA7SK`[1], HEXIM1 = ds$HEXIM1[1], DHFR3 = ds$DHFR3[1], DHFR12 = ds$DHFR12[1], meaning = ds$meaning[1] # 提取对应的标签信息
      )
    }))
    write_xlsx(gr_all, path = "output/growth_rate.xlsx")
  }
  # 对中值时间及生长速率进行线性拟合
  # tvfit <- lm(formula = time_odmid ~ growth_rate, data = gr_all)
  # gr_all$tmid_fitted <- fitted(tvfit)
  # r_squared <- summary(tvfit)$r.squared
  # print(summary(tvfit))
  # 对数拟合：t=(1/b)*ln(v/a*b)
  fit_line <- lm(formula = time_odmid ~ log(growth_rate), data = gr_all)
  gr_all$tmid_fitted <- fitted(fit_line)
  # 计算残差平方和 (RSS)
  rss <- sum((gr_all$time_odmid - gr_all$tmid_fitted)^2)
  # 计算总平方和 (TSS)
  tss <- sum((gr_all$time_odmid - mean(gr_all$time_odmid))^2)
  # 计算 R²
  r_squared <- 1 - rss / tss # 或等价于: r_squared <- sum((gr_all$tmid_fitted - y_mean)^2) / tss
  gr_all$class_color <- paste(gr_all$DHFR3, gr_all$DHFR12, sep = " ")
  gr_all$class_text <- paste(gr_all$meaning, gr_all$class_color, sep = " ")
  sample <- gr_all %>% # 分类取样
    group_by(class_text) %>% # 按class分组
    slice(1) %>% # 取每组第一行
    ungroup() # 取消分组
  # 生成高对比度颜色
  main_colors <- darken(c(
    `DHFR3_C_term DHFR12_N_term` = "#CABEEE",
    `DHFR3_C_term DHFR12_C_term` = "#ACE787",
    `DHFR3_N_term DHFR12_N_term` = "#B7FEF9",
    `DHFR3_N_term DHFR12_C_term` = "#F5C294",
    `NA NA` = "#F287F0"
  ), amount = 0.4, fix = TRUE)
  if (FALSE) {
    main_colors <- distinctColorPalette(length(unique(gr_all$class_color)), runTsne = TRUE)
    names(main_colors) <- unique(gr_all$class_color)
    dput(main_colors)
  }
  p_tv <- ggplot(
    data = gr_all,
    aes(x = growth_rate, y = time_odmid, color = class_color, shape = meaning)
  ) +
    geom_point() + # 实验数据点
    geom_label_repel(
      data = sample, aes(label = class_text),
      size = 2, alpha = 0.64, point.padding = 2, fill = NA,
    ) + # 图中添加数据标签
    geom_line(aes(x = growth_rate, y = tmid_fitted), # 拟合曲线
      color = "red", linetype = "dotdash", inherit.aes = FALSE
    ) +
    annotate("text",
      x = min(gr_all$growth_rate) + 0.15, y = max(gr_all$time_odmid) - 10,
      label = paste(
        "y=", round(fit_line$coefficients[1], 4), "+", round(fit_line$coefficients[2], 4), "*ln(x)", "\n",
        "R² =", round(r_squared, 3)
      ), hjust = 1
    ) +
    scale_color_manual(values = main_colors)
  ggsave("output/p_tv.png", plot = p_tv)
  # 类别与生长速率的关系
  p_vplasmid <- ggplot(
    data = gr_all,
    aes(x = growth_rate, y = plasmid, color = class_color, shape = meaning)
  ) +
    geom_point() + # 实验数据点
    geom_label_repel(
      data = sample, aes(label = class_text),
      size = 2, alpha = 0.64, point.padding = 2, fill = NA,
    ) + # 图中添加数据标签
    scale_color_manual(values = main_colors)
  ggsave("output/p_vp.png", plot = p_vplasmid)
  p_vp1 <- p_vplasmid + geom_boxplot(
    aes(group = class_text),
    fill = NA, alpha = 0.64,
    outlier.color = "red", outlier.shape = 4 # 离群值
  )
  ggsave("output/p_vp1.png", plot = p_vp1)
}
#----罗吉斯特模型拟合求生长曲线----
if (TRUE) {
  if (file.exists("output/growth_curve_logistic.xlsx")) {
    print("growth_curve_logistic.xlsx exists, loading it...")
    gr_lg <- read_excel("output/growth_curve_logistic.xlsx")
  } else {
    dir.create(path = "output", showWarnings = FALSE) # 创建输出目录
    gr_lg <- do.call("rbind", lapply(unique(od_long$well), function(x) {
      ds <- od_long[od_long$well == x, ]
      max_od <- max(ds$value)
      mid_od <- max_od - ((max_od - min(ds$value)) / 2)
      time_odmax <- ds$time[ds$value == max_od]
      time_odmid <- max(ds$time[ds$value <= mid_od & ds$time <= time_odmax[1]]) # mid_od不一定存在，且存在衰退期使odmax不止存在一个值，故取第一个
      time_odstart <- min(ds$time[ds$value > 5e-3]) # 取大于5e-3的最小时间点
      data_fitting <- data.frame(
        time_fitting = ds$time[ds$time <= 2 * time_odmid - time_odstart & ds$time >= time_odstart & ds$time >= 5e-3]
      ) # 从 ds 数据框中筛选出时间在 time_odmid 前后大于5e-3的数据，并将这些时间点存储到新的数据框 lmds 中，作为 t 列。
      data_fitting$value_fitting <- ds$value[ds$time %in% data_fitting$time_fitting]
      # data_fitting$value_fitting[data_fitting$value_fitting < 5e-3] <- 5e-3 # 将小于5e-3的值替换为5e-3
      model <- tryCatch(
        {
          drm(value_fitting ~ time_fitting, data = data_fitting, fct = L.5(fixed = c(NA, 5e-3, NA, NA, NA))) # L.5()函数用于定义五参数逻辑斯蒂模型
        },
        error = function(e) {
          message(paste("!!!Error in drm for experiment:", x, " - ", e$message))
          return(NULL)
        }
      )
      if (is.null(model)) {
        print(paste("Model fitting failed for well:", x))
        return(NULL) # 如果模型拟合失败，跳过当前数据集
      }
      p_values <- coef(summary(model))[, "p-value"]
      # 生成预测数据（100个点用于平滑曲线）
      pred_data <- data.frame(
        time_fitting = seq(min(data_fitting$time_fitting), max(data_fitting$time_fitting), length = 100)
      )

      # 计算预测值和置信区间（注意：predict.drc返回矩阵）
      pred <- predict(model,
        newdata = pred_data,
        interval = "confidence", # 均值置信区间
        level = 0.95
      ) # 95%置信水平
      # 将预测结果合并到数据框
      pred_data$pred <- pred[, "Prediction"]
      pred_data$lower <- pred[, "Lower"]
      pred_data$upper <- pred[, "Upper"]
      # 计算R²
      r2 <- 1 - (sum(residuals(model)^2) / (length(data_fitting$value_fitting) - 1 - 1)) /
        (sum((data_fitting$value_fitting - mean(data_fitting$value_fitting))^2) / (length(data_fitting$value_fitting) - 1)) # 计算调整后R²值
      # 绘制图形
      p <- tryCatch(
        {
          ggplot() +
            # 1. 原始数据点
            geom_point(
              data = data_fitting,
              aes(x = time_fitting, y = value_fitting),
              color = "black"
            ) +
            # 2. 拟合曲线
            geom_line(
              data = pred_data,
              aes(x = time_fitting, y = pred),
              color = "red",
              linewidth = 1.5,
              linetype = "dotdash"
            ) +
            # 3. 置信区间包络线
            geom_ribbon(
              data = pred_data,
              aes(x = time_fitting, ymin = lower, ymax = upper),
              fill = "blue",
              alpha = 0.2
            ) +
            # 4. 添加R²标签
            annotate(
              "text",
              x = min(data_fitting$time_fitting),
              y = max(data_fitting$value_fitting),
              label = sprintf("R² = %.4f", r2),
              hjust = 0,
              vjust = 1,
              size = 5,
              color = "darkred",
              fontface = "bold"
            ) +
            # 美化主题和标签
            labs(
              title = "Logistic Model - " %>% paste(x, ds$plasmid[1]),
              x = "time/h",
              y = "OD_600",
            )
        },
        error = function(e) {
          message(paste("!!!Error in printing figure:", x, " - ", e$message))
          return(NULL)
        }
      )
      print(paste(x, ":"))
      tryCatch(
        {
          ggsave(paste("output/logistic model_", x, ds$plasmid[1], ".png"), plot = p)
        },
        error = function(e) {
          message(paste("!!!Error in saving figure:", x, " - ", e$message))
        }
      )
      params <- coef(model) # 提取参数
      # plot(model, log = "y", main = paste("Logistic Model -", x), xlab = "Time/h", ylab = "OD_600")
      data.frame(
        well = x, plasmid = ds$plasmid[1], replicate = ds$replicate[1], meaning = ds$meaning[1],
        `RNA7SK` = ds$`RNA7SK`[1], HEXIM1 = ds$HEXIM1[1], DHFR3 = ds$DHFR3[1], DHFR12 = ds$DHFR12[1], # 提取对应的标签信息
        growth_rate = -1 * params[1], value_min = 5e-3, value_max = params[2], f = params[4], # 提取模型参数
        time_odmid = time_odmid, time_odmid_fitted = params[3], # 提取时间参数
        row.names = NULL
      )
    }))
    write_xlsx(gr_lg, path = "output/growth_curve_logistic.xlsx")
  }
}

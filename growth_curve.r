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
    cols = -c(
      "well", "plasmid", "replicate",
      "RNA7SK", "HEXIM1", "DHFR3", "DHFR12", "meaning"
    ),
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
if (TRUE) {
  # 生成数据
  od_long$plasmidid <- paste(od_long$plasmid, od_long$meaning)
  data_control <- od_long %>% filter(meaning %in% c("positive", "negative"))
  data_control <- subset(
    data_control,
    select = -c(RNA7SK, HEXIM1, DHFR3, DHFR12)
  )
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
    main_colors <- distinctColorPalette(
      length(unique(data_figure$plasmidid)),
      runTsne = TRUE
    )
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
        aes(
          x = time, y = value, group = well,
          color = plasmidid, alpha = replicate
        )
      ) +
        facet_grid(DHFR3 + DHFR12 ~ `RNA7SK` + `HEXIM1`) +
        geom_line() +
        scale_color_manual(values = all_colors) +
        theme(theme_classic())
    },
    error = function(e) {
      message(paste("!!!Error in p_1:", e$message))
    }
  )
  ggsave("output/01.png", plot = p_1, width = 12.36, height = 7.64)
  p_2 <- tryCatch(
    {
      p_1 + # 绘制分面对应的对照组曲线
        geom_line(
          data = data_control,
          aes(
            x = time, y = value, group = well,
            color = plasmidid, alpha = replicate
          ),
          inherit.aes = FALSE,
          linewidth = 0.1, linetype = "dashed"
        )
    },
    error = function(e) {
      message(paste("!!!Error in p_2:", e$message))
    }
  )
  ggsave("output/02.png", plot = p_2, width = 12.36, height = 7.64)
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
      od_max <- max(ds$value)
      od_mid <- od_max - ((od_max - min(ds$value)) / 2)
      time_odmax <- ds$time[ds$value == od_max]
      time_odmid <- max(
        ds$time[ds$value <= od_mid & ds$time <= time_odmax[1]]
      ) # od_mid不一定存在，且存在衰退期使odmax不止存在一个值，故取第一个
      lmds <- data.frame(
        t = ds$time[ds$time <= time_odmid + 4 & ds$time >= time_odmid - 4]
      ) # 从 ds 数据框中筛选出时间在 time_odmid 前后4h内的数据，并将这些时间点存储到新的数据框 lmds 中，作为 t 列。
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
              "y=", round(lmfit$coefficients[1], 4), "+",
              round(lmfit$coefficients[2], 4), "*x", "\n",
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
      }
      data.frame(
        well = x, plasmid = ds$plasmid[1], replicate = ds$replicate[1],
        growth_rate = lmfit$coefficients[2], time_odmid = time_odmid, od_max = od_max,
        `RNA7SK` = ds$`RNA7SK`[1], HEXIM1 = ds$HEXIM1[1],
        DHFR3 = ds$DHFR3[1], DHFR12 = ds$DHFR12[1], meaning = ds$meaning[1] # 提取对应的标签信息
      )
    }))
    write_xlsx(gr_all, path = "output/growth_rate.xlsx")
  }
  #----线性拟合结果绘图----
  # 图中分类标签取样
  gr_all$DHFR_position <- paste(gr_all$DHFR3, gr_all$DHFR12, sep = " ")
  gr_all$plasmid_class <- paste(gr_all$DHFR_position, gr_all$meaning, sep = " ")
  sample <- gr_all %>%
    group_by(plasmid_class) %>% # 按class分组
    slice(1) %>% # 取每组第一行
    ungroup() # 取消分组
  # 生成高对比度颜色
  main_colors <- darken(c(
    `DHFR3_C_term DHFR12_N_term double negative` = "#C8A2E3",
    `DHFR3_C_term DHFR12_C_term double negative` = "#FCF8F6",
    `DHFR3_N_term DHFR12_N_term double negative` = "#4FCAC1",
    `DHFR3_N_term DHFR12_C_term double negative` = "#82ABCC",
    `DHFR3_C_term DHFR12_N_term double positive` = "#A125B6",
    `DHFR3_N_term DHFR12_N_term double positive` = "#CC7D99",
    `DHFR3_C_term DHFR12_C_term double positive` = "#B4DE5E",
    `DHFR3_N_term DHFR12_C_term double positive` = "#86E29D",
    `DHFR3_C_term DHFR12_N_term negative` = "#E0CDA5",
    `DHFR3_C_term DHFR12_N_term positive` = "#F6AE71"
  ), amount = 0.4, fix = TRUE)
  # 自动生成颜色
  if (FALSE) {
    main_colors <- distinctColorPalette(
      length(unique(gr_all$plasmid_class)),
      runTsne = TRUE
    )
    names(main_colors) <- unique(gr_all$plasmid_class)
    dput(main_colors)
  }
  #----拟合参数再拟合----
  if (FALSE) {
    # 对中值时间及生长速率进行线性拟合
    if (FALSE) { # 线性拟合
      tvfit <- lm(formula = time_odmid ~ growth_rate, data = gr_all)
      gr_all$tmid_fitted <- fitted(tvfit)
      r_squared <- summary(tvfit)$r.squared
    }
    if (FALSE) { # 对数拟合：t=(1/b)*ln(v/a*b)
      tvfit <- lm(formula = time_odmid ~ log(growth_rate), data = gr_all)
      gr_all$tmid_fitted <- fitted(tvfit)
      r_squared <- summary(tvfit)$r.squared
    }
    p_tv <- ggplot(
      data = gr_all,
      aes(x = growth_rate, y = time_odmid, color = DHFR_position, shape = meaning)
    ) +
      geom_point() + # 实验数据点
      geom_label_repel(
        data = sample, aes(label = plasmid_class),
        size = 2, alpha = 0.64, point.padding = 2, fill = NA,
      ) + # 图中添加数据标签
      geom_line(aes(x = growth_rate, y = tmid_fitted), # 拟合曲线
        color = "red", linetype = "dotdash", inherit.aes = FALSE
      ) +
      annotate("text",
        x = min(gr_all$growth_rate) + 0.15, y = max(gr_all$time_odmid) - 10,
        label = paste(
          "y=", round(tvfit$coefficients[1], 4), "+",
          round(tvfit$coefficients[2], 4), "*x", "\n",
          "R² =", round(r_squared, 3)
        ), hjust = 1
      ) +
      scale_color_manual(values = main_colors)
    ggsave("output/p_tv.png", plot = p_tv, width = 12.36, height = 7.64)
    # 最大OD值与生长速率的关系
    odvfit <- lm(formula = od_max ~ growth_rate, data = gr_all)
    gr_all$odmax_fitted <- fitted(odvfit)
    r_squared <- summary(odvfit)$r.squared
    p_odv <- ggplot(
      data = gr_all,
      aes(x = growth_rate, y = od_max, color = DHFR_position, shape = meaning)
    ) +
      geom_point() + # 实验数据点
      geom_label_repel(
        data = sample, aes(label = plasmid_class),
        size = 2, alpha = 0.64, point.padding = 2, fill = NA,
      ) + # 图中添加数据标签
      geom_line(aes(x = growth_rate, y = odmax_fitted), # 拟合曲线
        color = "red", linetype = "dotdash", inherit.aes = FALSE
      ) +
      annotate("text",
        x = min(gr_all$growth_rate) + 0.05, y = max(gr_all$od_max),
        label = paste(
          "y=", round(tvfit$coefficients[1], 4), "+",
          round(tvfit$coefficients[2], 4), "*x", "\n",
          "R² =", round(r_squared, 3)
        ), hjust = 1
      ) +
      scale_color_manual(values = main_colors)
    ggsave("output/p_odv.png", plot = p_odv, width = 12.36, height = 7.64)
  }
  #----类别与生长速率的关系----
  p_vplasmid <- ggplot(
    data = gr_all,
    aes(x = plasmid, y = growth_rate, color = plasmid_class, shape = meaning)
  ) +
    facet_wrap(~DHFR_position, scales = "free_x", nrow = 1) + # 分组绘图
    geom_point() + # 实验数据点
    # geom_label_repel(
    #  data = sample, aes(label = plasmid_class),
    #  size = 2, alpha = 0.64, point.padding = 2, fill = NA,
    # ) + # 图中添加数据标签
    scale_color_manual(values = main_colors) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  ggsave("output/p_vp.png", plot = p_vplasmid, width = 12.36, height = 7.64)
  p_vp1 <- p_vplasmid + geom_boxplot(
    aes(group = plasmid_class),
    fill = NA, alpha = 0.64,
    outlier.color = "red", outlier.shape = 4 # 离群值
  )
  ggsave("output/p_vp1.png", plot = p_vp1, width = 12.36, height = 7.64)
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
      od_max <- max(ds$value)
      od_mid <- od_max - ((od_max - min(ds$value)) / 2)
      time_odmax <- ds$time[ds$value == od_max]
      time_odmid <- max(
        ds$time[ds$value <= od_mid & ds$time <= time_odmax[1]]
      ) # od_mid不一定存在，且存在衰退期使odmax不止存在一个值，故取第一个
      time_odstart <- min(ds$time[ds$value > 5e-3]) # 取大于5e-3的最小时间点
      data_fitting <- data.frame(
        time_fitting = ds$time[
          ds$time <= 2 * time_odmid - time_odstart & ds$time >= time_odstart &
            ds$time >= 5e-3
        ]
      ) # 从 ds 数据框中筛选出时间在 time_odmid 前后大于5e-3的数据，并将这些时间点存储到新的数据框 lmds 中，作为 t 列。
      data_fitting$value_fitting <- ds$value[ds$time %in% data_fitting$time_fitting]
      model <- tryCatch(
        {
          drm(value_fitting ~ time_fitting,
            data = data_fitting,
            fct = L.3(fixed = c(NA, NA, NA))
          ) # 逻辑斯蒂模型
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
      if (TRUE) { # 是否绘制详细拟合图
        p_values <- coef(summary(model))[, "p-value"]
        # 生成预测数据（100个点用于平滑曲线）
        pred_data <- data.frame(
          time_fitting = seq(
            min(data_fitting$time_fitting), max(data_fitting$time_fitting),
            length = 100
          )
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
        r2 <- 1 - (
          sum(residuals(model)^2) / (length(data_fitting$value_fitting) - 1 - 1)
        ) /
          (
            sum((data_fitting$value_fitting - mean(data_fitting$value_fitting))^2) / (length(data_fitting$value_fitting) - 1)
          ) # 计算调整后R²值
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
      }
      params <- coef(model) # 提取参数
      data.frame(
        well = x, plasmid = ds$plasmid[1], replicate = ds$replicate[1], meaning = ds$meaning[1],
        `RNA7SK` = ds$`RNA7SK`[1], HEXIM1 = ds$HEXIM1[1],
        DHFR3 = ds$DHFR3[1], DHFR12 = ds$DHFR12[1], # 提取对应的标签信息
        od_max = od_max, time_odmid = time_odmid,
        growth_rate = -1 * params[1], od_max_fitted = params[2], time_odmid_fitted = params[3], # 提取模型参数
        od_min_fitted = params[2] / (1 + exp(-1 * params[1] * params[3])),
        row.names = NULL
      )
    }))
    write_xlsx(gr_lg, path = "output/growth_curve_logistic.xlsx")
  }
  #----罗吉斯特模型拟合结果绘图----
  gr_lg$DHFR_position <- paste(gr_lg$DHFR3, gr_lg$DHFR12, sep = " ")
  gr_lg$plasmid_class <- paste(gr_lg$meaning, gr_lg$DHFR_position, sep = " ")
  sample <- gr_lg %>% # 分类取样
    group_by(plasmid_class) %>% # 按class分组
    slice(1) %>% # 取每组第一行
    ungroup() # 取消分组
  # 生成高对比度颜色
  main_colors <- darken(c(
    `double negative DHFR3_C_term DHFR12_N_term` = "#94B1DC",
    `double negative DHFR3_C_term DHFR12_C_term` = "#B1F1FA",
    `double negative DHFR3_N_term DHFR12_N_term` = "#82D45C",
    `double negative DHFR3_N_term DHFR12_C_term` = "#DED860",
    `double positive DHFR3_C_term DHFR12_N_term` = "#E573F8",
    `double positive DHFR3_N_term DHFR12_N_term` = "#C97A7A",
    `double positive DHFR3_C_term DHFR12_C_term` = "#AE88F3",
    `double positive DHFR3_N_term DHFR12_C_term` = "#B6C498",
    `negative DHFR3_C_term DHFR12_N_term` = "#8BEAB5",
    `positive DHFR3_C_term DHFR12_N_term` = "#F5D1F5"
  ), amount = 0.4, fix = TRUE)
  # 自动生成颜色
  if (FALSE) {
    main_colors <- distinctColorPalette(
      length(unique(gr_lg$plasmid_class)),
      runTsne = TRUE
    )
    names(main_colors) <- unique(gr_lg$plasmid_class)
    dput(main_colors)
  }
  #----拟合参数再拟合----
  if (FALSE) {
    # ----对中值时间及生长速率进行线性拟合----
    # gr_lg <- gr_lg[gr_lg$well != "K11", ] # 除去异常值
    tvfit <- lm(formula = time_odmid ~ growth_rate, data = gr_lg)
    gr_lg$tmid_fitted <- fitted(tvfit)
    r_squared <- summary(tvfit)$r.squared

    p_tv_lg <- ggplot(
      data = gr_lg,
      aes(x = growth_rate, y = time_odmid, color = DHFR_position, shape = meaning)
    ) +
      geom_point() + # 实验数据点
      geom_point(aes(x = growth_rate, y = time_odmid_fitted),
        alpha = 0.5
      ) + # 拟合数据点
      geom_label_repel(
        data = sample, aes(label = plasmid_class),
        size = 2, alpha = 0.64, point.padding = 2, fill = NA,
      ) + # 图中添加数据标签
      geom_line(aes(x = growth_rate, y = tmid_fitted), # 拟合曲线
        color = "red", linetype = "dotdash", inherit.aes = FALSE
      ) +
      annotate("text",
        x = min(gr_lg$growth_rate) + 0.3, y = max(gr_lg$time_odmid) - 10,
        label = paste(
          "y=", round(tvfit$coefficients[1], 4), "+",
          round(tvfit$coefficients[2], 4), "*x", "\n",
          "R² =", round(r_squared, 3)
        ), hjust = 1
      ) +
      scale_color_manual(values = main_colors)
    ggsave("output/p_tv_lg.png", plot = p_tv_lg, width = 12.36, height = 7.64)
    # ----最大OD值与生长速率进行线性拟合----
    odvfit <- lm(formula = od_max ~ growth_rate, data = gr_lg)
    gr_lg$odmax_fitted <- fitted(odvfit)
    r_squared <- summary(odvfit)$r.squared
    p_odv_lg <- ggplot(
      data = gr_lg,
      aes(x = growth_rate, y = od_max, color = DHFR_position, shape = meaning)
    ) +
      geom_point() + # 实验数据点
      geom_point(aes(x = growth_rate, y = od_max_fitted), alpha = 0.5) + # 拟合数据点
      geom_label_repel(
        data = sample, aes(label = plasmid_class),
        size = 2, alpha = 0.64, point.padding = 2, fill = NA,
      ) + # 图中添加数据标签
      geom_line(aes(x = growth_rate, y = odmax_fitted), # 拟合曲线
        color = "red", linetype = "dotdash", inherit.aes = FALSE
      ) +
      annotate("text",
        x = min(gr_lg$growth_rate) + 0.2, y = max(gr_lg$od_max),
        label = paste(
          "y=", round(tvfit$coefficients[1], 4), "+",
          round(tvfit$coefficients[2], 4), "*x", "\n",
          "R² =", round(r_squared, 3)
        ), hjust = 1
      ) +
      scale_color_manual(values = main_colors)
    ggsave("output/p_odv_lg.png", plot = p_odv_lg, width = 12.36, height = 7.64)
  }
  # ----类别与生长速率的关系----
  p_vplasmid_lg <- ggplot(
    data = gr_lg,
    aes(
      x = plasmid, y = growth_rate,
      color = plasmid_class, shape = meaning,
      alpha = replicate
    )
  ) +
    facet_wrap(~DHFR_position, scales = "free_x", nrow = 1) + # 分组绘图
    geom_point() + # 实验数据点
    # geom_label_repel(
    #  data = sample, aes(label = plasmid_class),
    #  size = 2, alpha = 0.64, point.padding = 2, fill = NA,
    # ) + # 图中添加数据标签
    scale_color_manual(values = main_colors) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  ggsave("output/p_vp_lg.png", plot = p_vplasmid_lg, width = 12.36, height = 7.64)
  p_vp1_lg <- p_vplasmid_lg + geom_boxplot(
    aes(group = plasmid_class),
    fill = NA, alpha = 0.64,
    outlier.color = "red", outlier.shape = 4 # 离群值
  )
  ggsave("output/p_vp1_lg.png", plot = p_vp1_lg, width = 12.36, height = 7.64)
  # ----类别与od_max的关系----
  p_odmaxplasmid_lg <- ggplot(
    data = gr_lg,
    aes(
      x = plasmid, y = od_max,
      color = plasmid_class, shape = meaning,
      alpha = replicate
    )
  ) +
    facet_wrap(~DHFR_position, scales = "free_x", nrow = 1) + # 分组绘图
    geom_point() + # 实验数据点
    # geom_label_repel(
    #  data = sample, aes(label = plasmid_class),
    #  size = 2, alpha = 0.64, point.padding = 2, fill = NA,
    # ) + # 图中添加数据标签
    scale_color_manual(values = main_colors) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  ggsave("output/p_odp_lg.png", plot = p_odmaxplasmid_lg, width = 12.36, height = 7.64)
  p_odp1_lg <- p_odmaxplasmid_lg + geom_boxplot(
    aes(group = plasmid_class),
    fill = NA, alpha = 0.64,
    outlier.color = "red", outlier.shape = 4 # 离群值
  )
  ggsave("output/p_odp1_lg.png", plot = p_odp1_lg, width = 12.36, height = 7.64)
}

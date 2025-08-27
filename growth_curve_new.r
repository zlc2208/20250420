library(readxl)
library(writexl)
library(crayon)
library(dplyr)
library(tidyr)
library(drc)
library(tibble)
library(ggplot2)

theme_set(theme_bw())
date <- as.character("20250420")

#----读取数据并处理为长格式----
if (file.exists(paste(
  "output_", date, "/od_long_", date, ".xlsx",
  sep = ""
))) {
  cat(crayon::green("od_long.xlsx exists, loading it...\n"))
  od_long <- read_excel(paste(
    "output_", date, "/od_long_", date, ".xlsx",
    sep = ""
  )) # 读取数据
} else {
  cat(crayon::red("od_long.xlsx does not exist, creating it...\n"))
  ## 给OD值加上标签并转化为长格式:od_long
  if (TRUE) {
    OD_xlsx_longer <- function(file_data = "data_date.xlsx", file_label = "label_date.xlsx") {
      data_od <- read_excel(file_data) %>%
        column_to_rownames(var = "well")
      label <- read_excel(file_label, sheet = "label") # 读取标签
      # OD值减去blank均值
      OD_0 <- data_od %>%
        filter(row.names(data_od) %in% label$well[label$plasmid == "blank"]) %>%
        summarise(across(where(is.numeric), mean)) %>%
        as.numeric()
      real_od <- t(apply(data_od, 1, function(i) {
        i - OD_0
      })) %>%
        data.frame()
      colnames(real_od) <- colnames(data_od) # 恢复列名
      real_od <- rownames_to_column(real_od, var = "well")
      od_melt <- merge(real_od, label, by = "well")
      od_melt <- od_melt[!(od_melt$plasmid == "blank"), ]
      num_cols <- names(od_melt)[suppressWarnings(!is.na(as.numeric(names(od_melt))))]
      od_long <- pivot_longer(
        od_melt,
        cols = all_of(num_cols),
        names_to = "time",
        values_to = "value",
        names_transform = list(time = as.numeric),
        values_transform = list(value = as.numeric)
      )
      od_long$time <- sapply(od_long$time, function(x) {
        x / 6
      })
      od_long$value[od_long$value < 0] <- 0
      dir.create(paste(
        "output_", date,
        sep = ""
      ), showWarnings = FALSE, recursive = TRUE)
      write_xlsx(od_long,
        path = paste("output_", date, "/od_long_", date, ".xlsx", sep = "")
      )
      return(od_long)
    }
    save(OD_xlsx_longer, file = "d:/app/r/function/OD_xlsx_longer.RData")
  } else {
    load("d:/app/r/function/OD_xlsx_longer.RData")
  }
  od_long <- OD_xlsx_longer(
    file_data = paste("data_", date, ".xlsx", sep = ""),
    file_label = paste("label_", date, ".xlsx", sep = "")
  )
}

#----罗吉斯特模型拟合求生长曲线----
if (TRUE) {
  if (file.exists(paste(
    "output_", date, "/growth_curve_logistic_", date, ".xlsx",
    sep = ""
  ))) {
    cat(crayon::green("growth_curve_logistic.xlsx exists, loading it...\n"))
    gr_lg <- read_excel(paste(
      "output_", date, "/growth_curve_logistic_", date, ".xlsx",
      sep = ""
    ))
  } else {
    cat(crayon::red("growth_curve_logistic.xlsx does not exist, creating it...\n"))
    if (TRUE) {
      OD_longer_grlg <- function(od_long = od_long) {
        gr_lg <- do.call("rbind", lapply(unique(od_long$well), function(x) {
          gr_lg_sub <- od_long[od_long$well == x, ]
          od_max <- max(gr_lg_sub$value)
          od_mid <- od_max - ((od_max - min(gr_lg_sub$value)) / 2)
          time_odmax <- gr_lg_sub$time[gr_lg_sub$value == od_max]
          time_odmid <- max(
            gr_lg_sub$time[gr_lg_sub$value <= od_mid & gr_lg_sub$time <= time_odmax[1]]
          ) # od_mid不一定存在，且存在衰退期使odmax不止存在一个值，故取第一个
          time_start <- min(gr_lg_sub$time[gr_lg_sub$value >= 5e-3]) # 取大于5e-3的最小时间点
          data_fitting <- data.frame(
            time_fitting = gr_lg_sub$time[gr_lg_sub$time <= 2 * time_odmid - time_start & gr_lg_sub$time >= time_start]
          ) # 从 gr_lg_sub 数据框中筛选出时间在 time_odmid 前后大于5e-3的数据，并将这些时间点存储到新的数据框 gr_lg_sub 中，作为 t 列。
          data_fitting$value_fitting <- gr_lg_sub$value[
            gr_lg_sub$time %in% data_fitting$time_fitting
          ]
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
            cat(crayon::red(paste("Model fitting failed for well:", x, "\n")))
            # 如果模型拟合失败，跳过当前数据集
          } else {
            params <- coef(model) # 提取参数
            # 是否绘制详细拟合图
            if (TRUE) {
              # calculate more data
              time_fitted <- data.frame(time_fitting = gr_lg_sub$time)
              pred <- predict(model, newdata = time_fitted, interval = "confidence", level = 0.95)
              data_fitted <- cbind(time_fitted, pred)
              colnames(data_fitted) <- c("time_fitted", "value_fitted", "lwr", "upr")
              # 绘制图形
              p <- ggplot(data_fitted, aes(x = time_fitted, y = value_fitted)) +
                geom_point( # original data point
                  data = gr_lg_sub,
                  aes(x = time, y = value, color = "original data"),
                  inherit.aes = FALSE
                ) + # fit data
                geom_line(
                  aes(color = "fit line"),
                  linewidth = 1.5, linetype = "dotdash"
                ) + # confidence interval
                geom_ribbon(
                  aes(ymin = lwr, ymax = upr),
                  fill = "blue", alpha = 0.2
                ) + # fit interval
                geom_vline(
                  xintercept = c(min(data_fitting$time_fitting), max(data_fitting$time_fitting)),
                  color = "red", linetype = "dotted", linewidth = 1.2
                ) +
                annotate("text",
                  label = "Fit Interval", x = max(data_fitting$time_fitting), y = min(data_fitting$value_fitting),
                  hjust = 1.1, vjust = -0.1, color = "red"
                ) +
                scale_color_manual(values = c("original data" = "black", "fit line" = "red")) +
                labs(
                  title = paste("Logistic Model - ", x, gr_lg_sub$plasmidid[1]),
                  x = "time/h", y = expression(OD[600]), color = "logistic model fit"
                )
              # add regulated r2
              n_fitting <- length(data_fitting$value_fitting)
              r2 <- 1 - (sum(residuals(model)^2) / (n_fitting - length(params) - 1)) /
                (sum((data_fitting$value_fitting - mean(data_fitting$value_fitting))^2) / (n_fitting - 1))
              p <- p + annotate("text",
                label = bquote(OD[600] == frac(
                  .(round(params[2], 3)),
                  1 + e^{
                    .(round(params[1], 3)) * (time - .(round(params[3], 3)))
                  }
                )),
                x = min(data_fitting$time_fitting), y = max(data_fitting$value_fitting),
                hjust = -0.1, vjust = 1, size = 5, color = "black", fontface = "bold"
              ) + annotate("text",
                label = bquote(r^2 == .(round(r2, 4))),
                x = min(data_fitting$time_fitting), y = max(data_fitting$value_fitting),
                hjust = -0.1, vjust = 4, size = 5, color = "darkred", fontface = "bold"
              )
              print(paste("logistic model fitting for well", x, "has done."))
              ggsave(
                paste(
                  "output_", date, "/logistic model_",
                  x, "_", gr_lg_sub$plasmidid[1], "_", date, ".png",
                  sep = ""
                ),
                plot = p, width = 8, height = 6
              )
            }
            data.frame( # get fit parameters
              od_max = od_max, time_odmid = time_odmid, adjusted_R2 = r2,
              growth_rate = -1 * params[1], od_max_fitted = params[2],
              time_odmid_fitted = params[3],
              od_min_fitted = params[2] / (1 + exp(-1 * params[1] * params[3]))
            ) %>%
              cbind(subset(gr_lg_sub[1, ], select = -c(time, value))) # get original label
          }
        }))
        write_xlsx(gr_lg, path = paste(
          "output_", date, "/growth_curve_logistic_", date, ".xlsx",
          sep = ""
        ))
        return(gr_lg)
      }
      save(OD_longer_grlg, file = "d:/app/r/function/OD_longer_grlg.RData")
    } else {
      load("d:/app/r/function/OD_longer_grlg.RData")
    }
    od_long$plasmidid <- paste(od_long$plasmid, od_long$meaning)
    gr_lg <- OD_longer_grlg(od_long = od_long)
  }
}

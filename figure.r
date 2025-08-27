library(readxl)
library(crayon)
library(ggplot2)
library(dplyr)
library(tidyr)
library(writexl)
library(randomcoloR)
library(gridExtra)
library(ggrepel)
library(colorspace)
library(drc)
library(tibble)
library(ggsignif)
library(rlang)
library(pheatmap)
library(RColorBrewer)
library(effsize)

theme_set(theme_bw())
date <- as.character("20250420")

#----loading growth curve data as od_long----
if (TRUE) {
    if (file.exists(paste(
        "output_", date, "/od_long_", date, ".xlsx",
        sep = ""
    ))) {
        cat(crayon::green("od_long.xlsx exists, loading it...\n"))
        od_long <- read_excel(paste("output_", date, "/od_long_", date, ".xlsx", sep = "")) %>%
            mutate(meaning = factor(meaning, levels = c("double positive", "double negative", "positive", "negative")))
    } else {
        cat(crayon::red("od_long.xlsx does not exist\n"))
    }
}
#----loading logistic fit data as gr_lg----
if (TRUE) {
    if (file.exists(paste(
        "output_", date, "/growth_curve_logistic_", date, ".xlsx",
        sep = ""
    ))) {
        cat(crayon::green("growth_curve_logistic.xlsx exists, loading it...\n"))
        gr_lg <- read_excel(paste("output_", date, "/growth_curve_logistic_", date, ".xlsx", sep = "")) %>%
            mutate(meaning = factor(meaning, levels = c("double positive", "double negative", "positive", "negative")))
    } else {
        cat(crayon::red("growth_curve_logistic.xlsx does not exist\n"))
    }
}

#----growth curve
if (FALSE) {
    # 生成数据
    od_long$plasmidid <- paste(od_long$plasmid, od_long$meaning)
    data_control <- od_long %>% filter(meaning %in% c("positive", "negative"))
    data_control <- subset(
        data_control,
        select = -c(RNA, Protein, DHFR_layout)
    )
    data_figure <- od_long %>% filter(!meaning %in% c("positive", "negative"))
    # ----绘图----
    # 步骤1：生成主图的自动颜色
    if (FALSE) {
        color_figure <- distinctColorPalette(
            length(unique(data_figure$plasmidid)),
            runTsne = TRUE
        )
        names(color_figure) <- unique(data_figure$plasmidid)
        cat(crayon::yellow("color_figure generated successfully.\n"))
        dput(color_figure)
    } else {
        color_figure <- darken(c(
            `pTY126-1 double positive` = "#91C6F4",
            `pTY126-3 double positive` = "#4DCEAD",
            `pTY126-2 double positive` = "#86E9F4",
            `pTY126-4 double positive` = "#F5ACE4",
            `pTY127-1 double positive` = "#A8DFAA",
            `pTY127-3 double positive` = "#8070D7",
            `pTY127-2 double positive` = "#CBCCAB",
            `pTY127-4 double positive` = "#8CFBB2",
            `pTY128-1 double positive` = "#BAAAE6",
            `pTY128-3 double positive` = "#D2B55B",
            `pTY128-2 double positive` = "#EF7190",
            `pTY128-4 double positive` = "#90BAB0",
            `pTY129-1 double positive` = "#D2EED9",
            `pTY129-3 double positive` = "#B09282",
            `pTY129-2 double positive` = "#E3F560",
            `pTY129-4 double positive` = "#9F20C8",
            `pTY102-1 double negative` = "#9D2B99",
            `pTY102-3 double negative` = "#6CCA54",
            `pTY102-2 double negative` = "#D7ECA0",
            `pTY102-4 double negative` = "#DA8657",
            `pTY120-1 double negative` = "#558CD0",
            `pTY121-5 double negative` = "#CBD4E2",
            `pTY121-1 double negative` = "#CE82E2",
            `pTY104-1 double negative` = "#FED2E3"
        ), amount = 0.382, fix = TRUE)
    }
    # 步骤2：手动设置对照组颜色
    if (FALSE) {
        color_control <- distinctColorPalette(
            length(unique(data_control$plasmidid)),
            runTsne = TRUE
        )
        names(color_control) <- unique(data_control$plasmidid)
        cat(crayon::yellow("color_control generated successfully.\n"))
        dput(color_control)
    } else {
        color_control <- c(
            `pGJJ588 positive` = "pink",
            `pGJJ336 positive` = "pink",
            `pGJJ584 negative` = "lightblue"
        )
    }
    # 步骤3：合并颜色，确保对照组颜色覆盖同名类别
    color_all <- color_figure
    color_all[names(color_control)] <- color_control
    p_1 <- ggplot(
        data = data_figure,
        aes(
            x = time, y = value, group = well,
            color = plasmidid, alpha = replicate
        )
    ) +
        facet_grid(DHFR_layout ~ meaning) +
        geom_line() +
        scale_color_manual(values = color_all) +
        labs(x = "time/h", y = expression(OD[600]))
    ggsave(paste(
        "figure_", date, "/growth_curve_", date, ".png",
        sep = ""
    ), plot = p_1, width = 11, height = 7.64)
    cat(crayon::green("Plot p_1 saved successfully.\n"))
    p_2 <- ggplot(
        data = data_figure,
        aes(
            x = time, y = value, group = well,
            color = plasmidid, alpha = replicate
        )
    ) +
        facet_grid(DHFR_layout ~ Protein + RNA) +
        geom_line() +
        scale_color_manual(values = color_all) +
        labs(x = "time/h", y = expression(OD[600])) + # 绘制分面对应的对照组曲线
        geom_line(
            data = data_control,
            aes(
                x = time, y = value, group = well,
                color = plasmidid, alpha = replicate
            ),
            inherit.aes = FALSE,
            linetype = "dashed"
        )
    ggsave(paste(
        "figure_", date, "/growth_curve_control_", date, ".png",
        sep = ""
    ), plot = p_2, width = 12.36, height = 7.64)
    cat(crayon::green("Plot p_2 saved successfully.\n"))
}

#----boxplot----
if (TRUE) {
    gr_lg$plasmid_type <- paste(gr_lg$meaning, gr_lg$DHFR_layout, sep = " ")
    gr_lg_figure <- subset(gr_lg, subset = !(gr_lg$meaning %in% c("positive", "negative")))
    # 自动生成颜色
    if (FALSE) {
        color_figure <- distinctColorPalette(
            length(unique(gr_lg_figure$plasmid_type)),
            runTsne = TRUE
        )
        names(color_figure) <- unique(gr_lg_figure$plasmid_type)
        cat(crayon::yellow("color_figure generated successfully.\n"))
        dput(color_figure)
    } else {
        color_figure <- darken(c(
            `double negative DHFR_12N3C` = "#88633B", `double negative DHFR_12C3C` = "#F2FD95",
            `double negative DHFR_12N3N` = "#14BFA6", `double negative DHFR_12C3N` = "#9BDCFD",
            `double positive DHFR_12N3C` = "#DCE8D0", `double positive DHFR_12N3N` = "#7AD56D",
            `double positive DHFR_12C3C` = "#E681C3", `double positive DHFR_12C3N` = "#D2C9FB"
        ), amount = 0.382, fix = TRUE)
    }
    # ----类别与od_max的关系----
    odmax_plasmid <- ggplot(
        data = gr_lg_figure,
        aes(x = plasmid, y = od_max, color = plasmid_type, shape = meaning, alpha = replicate)
    ) +
        facet_wrap(~DHFR_layout, scales = "free_x", nrow = 1) + # 分组绘图
        geom_point(size = 2) +
        geom_boxplot(
            aes(group = plasmid_type),
            alpha = 0.64, width = .8,
            outlier.color = "red", outlier.shape = 4 # 离群值
        ) +
        stat_boxplot(
            geom = "errorbar", aes(group = plasmid_type),
            alpha = .64, width = .43
        ) +
        scale_color_manual(values = color_figure) +
        labs(y = expression(OD[600])) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    ggsave(
        paste("figure_", date, "/odmax_plasmid_", date, ".png", sep = ""),
        plot = odmax_plasmid, width = 12.36, height = 7.64
    )
    cat(crayon::green("Plot odmax_plasmid saved successfully.\n"))
    # add p_value of t_test
    # test_list <- lapply(unique(gr_lg_figure$DHFR_layout), function(x) {
    #  plasmid_type <- unique(gr_lg_figure$plasmid_type[gr_lg_figure$DHFR_layout == x])
    #  combn(plasmid_type, 2, simplify = FALSE) %>%
    #    unlist(recursive = FALSE) %>%
    #    as.vector()
    # })
    odmax_plasmid_t <- ggplot(gr_lg_figure, aes(x = meaning, y = od_max, color = plasmid_type)) +
        facet_wrap(~DHFR_layout, nrow = 1) +
        geom_boxplot(outlier.color = "red", outlier.shape = 4, width = .8) +
        stat_boxplot(geom = "errorbar", width = .43) +
        geom_jitter(aes(shape = meaning), width = .3, size = 1.5) +
        scale_color_manual(values = color_figure) +
        labs(title = expression(
            "T-test: " * H[0] * ":" * OD600[max["double positive"]] * "=" * OD600[max["double negative"]]
        ), y = expression(OD[600])) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    # for (i in test_list) {
    #  data_i <- gr_lg_figure[gr_lg_figure$plasmid_type %in% i, ]
    #  odmax_plasmid_t <- odmax_plasmid_t +
    #    geom_signif(
    #      data = data_i,
    #      comparisons = list(i),
    #      test = "t.test", color = "darkred",
    #      map_signif_level = TRUE # 自动转换为星号 (*, **, ***)
    #    )
    # }
    for (i in unique(gr_lg_figure$DHFR_layout)) {
        data_i <- gr_lg_figure[gr_lg_figure$DHFR_layout %in% i, ]
        odmax_plasmid_t <- odmax_plasmid_t +
            geom_signif(
                data = data_i, test = "t.test", color = "darkred",
                comparisons = list(c("double negative", "double positive")),
                map_signif_level = TRUE,
                y_position = max(data_i$od_max)
            )
    }
    ggsave(paste("figure_", date, "/odmax_plasmid_t_", date, ".png", sep = ""),
        plot = odmax_plasmid_t, width = 12, height = 5
    )
    cat(crayon::green("Plot odmax_plasmid_t saved successfully.\n"))

    # ----类别与生长速率的关系----
    growth_rate_plasmid <- ggplot(
        data = gr_lg_figure,
        aes(x = plasmid, y = growth_rate, color = plasmid_type, shape = meaning, alpha = replicate)
    ) +
        facet_wrap(~DHFR_layout, scales = "free_x", nrow = 1) + # 分组绘图
        geom_point(size = 2) +
        geom_boxplot(
            aes(group = plasmid_type),
            alpha = 0.64, width = .8,
            outlier.color = "red", outlier.shape = 4 # 离群值
        ) +
        stat_boxplot(
            geom = "errorbar", aes(group = plasmid_type),
            alpha = .64, width = .43
        ) +
        scale_color_manual(values = color_figure) +
        labs(y = expression("Growth_rate/(" * OD[600] * "/h)")) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    ggsave(paste("figure_", date, "/growth_rate_plasmid_", date, ".png", sep = ""),
        plot = growth_rate_plasmid, width = 12.36, height = 7.64
    )
    cat(crayon::green("Plot growth_rate_plasmid saved successfully.\n"))
    # add p_value of t_test
    # test_list <- lapply(unique(gr_lg_figure$DHFR_layout), function(x) {
    # plasmid_type <- unique(gr_lg_figure$plasmid_type[gr_lg_figure$DHFR_layout == x])
    #  combn(plasmid_type, 2, simplify = FALSE) %>%
    #    unlist(recursive = FALSE) %>%
    #    as.vector()
    # })
    growth_rate_plasmid_t <- ggplot(gr_lg_figure, aes(x = meaning, y = growth_rate, color = plasmid_type)) +
        facet_wrap(~DHFR_layout, nrow = 1) +
        geom_boxplot(outlier.color = "red", outlier.shape = 4, width = .8) +
        stat_boxplot(geom = "errorbar", width = .43) +
        geom_jitter(aes(size = meaning), width = .3, size = 1.5) +
        scale_color_manual(values = color_figure) +
        labs(title = expression(
            "T-test: " * H[0] * ":" * Growth_rate["double positive"] * "=" * Growth_rate["double negative"]
        ), y = expression("Growth_rate/(" * OD[600] * "/h)")) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    # for (i in test_list) {
    #  data_i <- gr_lg_figure[gr_lg_figure$plasmid_type %in% i, ]
    #  growth_rate_plasmid_t <- growth_rate_plasmid_t +
    #    geom_signif(
    #      data = data_i,
    #      comparisons = list(i),
    #      test = "t.test", color = "darkred",
    #      map_signif_level = TRUE # 自动转换为星号 (*0.05, **0.01, ***0.001)
    #    )
    # }
    for (i in unique(gr_lg_figure$DHFR_layout)) {
        data_i <- gr_lg_figure[gr_lg_figure$DHFR_layout %in% i, ]
        growth_rate_plasmid_t <- growth_rate_plasmid_t +
            geom_signif(
                data = data_i, test = "t.test", color = "darkred",
                comparisons = list(c("double negative", "double positive")),
                map_signif_level = TRUE,
                y_position = max(data_i$growth_rate)
            )
    }
    ggsave(paste("figure_", date, "/growth_rate_plasmid_t_", date, ".png", sep = ""),
        plot = growth_rate_plasmid_t, width = 12, height = 5
    )
    cat(crayon::green("Plot growth_rate_plasmid_t saved successfully.\n"))
    #----bar----
    if (TRUE) {
        # odmax_bar
        odmax_bar <- ggplot(gr_lg_figure, aes(x = meaning, y = od_max, fill = plasmid_type)) +
            facet_wrap(~DHFR_layout, nrow = 1) +
            geom_jitter(aes(color = plasmid_type), fill = "white", width = 0.2, pch = 1) +
            scale_color_manual(values = color_figure) +
            stat_summary(fun = mean, geom = "bar", width = 0.6, alpha = 0.64) + # 绘制条状图
            stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # 添加误差线 (标准误)
            scale_fill_manual(values = color_figure) +
            labs(
                y = expression(OD[600]),
                title = expression(
                    "T-test: " * H[0] * ":" * OD600[max["double positive"]] * "=" * OD600[max["double negative"]]
                )
            ) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        for (i in unique(gr_lg_figure$DHFR_layout)) {
            odmax_bar <- odmax_bar +
                geom_signif(
                    data = gr_lg_figure[gr_lg_figure$DHFR_layout == i, ],
                    test = "t.test", color = "darkred",
                    comparisons = list(c("double positive", "double negative")),
                    map_signif_level = TRUE, step_increase = 0.1
                )
        }
        ggsave(paste("figure_", date, "/odmax_bar_", date, ".png", sep = ""),
            plot = odmax_bar, width = 10, height = 5
        )
        cat(crayon::green("Plot odmax_single_bar saved successfully.\n"))
        # growth_rate bar
        growth_rate_bar <- ggplot(gr_lg_figure, aes(x = meaning, y = growth_rate, fill = plasmid_type)) +
            facet_wrap(~DHFR_layout, nrow = 1) +
            geom_jitter(aes(color = plasmid_type), fill = "white", width = 0.2, pch = 1) +
            scale_color_manual(values = color_figure) +
            stat_summary(fun = mean, geom = "bar", width = 0.6, alpha = 0.64) + # 绘制条状图
            stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + # 添加误差线 (标准误)
            scale_fill_manual(values = color_figure) +
            labs(
                y = expression("Growth rate/(" * OD[600] * "/h)"),
                title = expression(
                    "T-test: " * H[0] * ":" * Growth_rate["double positive"] * "=" * Growth_rate["double negative"]
                )
            ) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        for (i in unique(gr_lg_figure$DHFR_layout)) {
            growth_rate_bar <- growth_rate_bar +
                geom_signif(
                    data = gr_lg_figure[gr_lg_figure$DHFR_layout == i, ],
                    test = "t.test", color = "darkred",
                    comparisons = list(c("double positive", "double negative")),
                    map_signif_level = TRUE, step_increase = 0.1
                )
        }
        ggsave(paste("figure_", date, "/growth_rate_bar_", date, ".png", sep = ""),
            plot = growth_rate_bar, width = 10, height = 5
        )
        cat(crayon::green("Plot growth_rate_single_bar saved successfully.\n"))
    }
}

#----heatmap----
if (FALSE) {
    heatmap_MaV <- function(data_hp = data_hp, variable = variable, rname = rname, cname = cname, type = type) {
        name_variable <- as_name(ensym(variable))
        rname_row <- as_name(ensym(rname))
        cname_col <- as_name(ensym(cname))
        name_type <- as_name(ensym(type))
        data_MaV <- do.call("rbind", lapply(unique(data_hp[[name_type]]), function(x) {
            data_hp_sub <- data_hp[data_hp[[name_type]] == x, ]
            data.frame(
                rname = unique(data_hp_sub[[rname_row]]),
                cname = unique(data_hp_sub[[cname_col]]),
                value_mean = mean(data_hp_sub[[name_variable]]),
                value_sd = sd(data_hp_sub[[name_variable]])
            )
        })) %>%
            mutate(value_CV = value_sd / value_mean * 100) %>%
            arrange(rname, cname)
        # 平均值矩阵
        mean_mat <- pivot_wider(subset(data_MaV, select = c(rname, cname, value_mean)),
            names_from = cname, values_from = value_mean
        ) %>%
            column_to_rownames("rname") %>%
            as.matrix()
        print(paste("Mean Matrix_", as.character(name_variable), ":", seq = ""))
        print(mean_mat)
        pheatmap(mean_mat,
            main = paste("Mean of ", as.character(name_variable), seq = ""),
            cluster_rows = FALSE, cluster_cols = FALSE,
            color = colorRampPalette(brewer.pal(9, "Blues")[1:7])(255),
            angle_col = 45, #  热图列名角度
            display_numbers = TRUE, #  矩阵的数值是否显示在热图上
            number_format = "%.4f",
            filename = paste("figure_", date, "/heatmap_mean_",
                as.character(name_variable), "_", date, ".png",
                sep = ""
            ), width = 4, height = 8
        )
        cat(green(paste("Mean Heatmap_", as.character(name_variable), " finish.\n")))
        # 变异系数矩阵（小10%中100%大）
        CV_mat <- pivot_wider(subset(data_MaV, select = c(rname, cname, value_CV)),
            names_from = cname,
            values_from = value_CV,
            values_fn = mean
        ) %>%
            column_to_rownames("rname") %>%
            as.matrix()
        print(paste("CV Matrix_", as.character(name_variable), ":", seq = ""))
        print(CV_mat)
        breaks <- c(seq(0, 34, length.out = 255), 35, seq(36, ifelse(max(data_MaV$value_CV) > 100, max(data_MaV$value_CV), 100), length.out = 255))
        colors <- c(colorRampPalette(rev(brewer.pal(9, "Blues")[1:6]))(255), "white", colorRampPalette(brewer.pal(9, "Reds")[1:7])(255))
        pheatmap(CV_mat,
            main = paste("CV% of ", as.character(name_variable), seq = ""),
            cluster_rows = FALSE, cluster_cols = FALSE,
            breaks = breaks,
            color = colors,
            angle_col = 45, #  热图列名角度
            display_numbers = TRUE, #  矩阵的数值是否显示在热图上
            number_format = "%.2f",
            filename = paste("figure_", date,
                "/heatmap_CV_", as.character(name_variable), "_", date, ".png",
                sep = ""
            ), width = 4, height = 8
        )
        cat(green(paste("CV Heatmap_", as.character(name_variable), " finish.\n")))
    }
    save(heatmap_MaV, file = "d:/app/r/function/heatmap_MaV.RDate")
} else {
    load("d:/app/r/function/heatmap_MaV.RDate")
}
# paint the heatmap
if (FALSE) {
    gr_lg_hp <- gr_lg %>% filter(!meaning %in% c("positive", "negative"))
    gr_lg_hp$plasmid_type <- paste(gr_lg_hp$meaning, gr_lg_hp$DHFR_layout)
    heatmap_MaV(data_hp = gr_lg_hp, variable = od_max, rname = DHFR_layout, cname = meaning, type = plasmid_type)
    heatmap_MaV(data_hp = gr_lg_hp, variable = growth_rate, rname = DHFR_layout, cname = meaning, type = plasmid_type)
    heatmap_MaV(data_hp = gr_lg_hp, variable = adjusted_R2, rname = DHFR_layout, cname = meaning, type = plasmid_type)
}

#----t_test matrix----
if (FALSE) {
    gr_lg$plasmid_type <- paste(gr_lg$meaning, gr_lg$DHFR_layout)
    data_test <- subset(gr_lg,
        subset = !(gr_lg$meaning %in% c("positive", "negative")),
        select = c(plasmid_type, od_max, growth_rate, meaning, DHFR_layout)
    )
    types <- unique(data_test$plasmid_type)
    diff_pairs <- combn(types, 2, simplify = FALSE)
    diff_pairs_rev <- lapply(diff_pairs, function(pair) rev(pair))
    self_pairs <- lapply(types, function(x) c(x, x))
    all_pairs <- c(diff_pairs, diff_pairs_rev, self_pairs)
    #---- odmax----
    t_odmax <- do.call(rbind, lapply(all_pairs, function(pair) {
        test <- t.test(
            data_test$od_max[data_test$plasmid_type == pair[1]],
            data_test$od_max[data_test$plasmid_type == pair[2]],
            alternative = "greater"
        )
        cd <- cohen.d(
            data_test$od_max[data_test$plasmid_type == pair[1]],
            data_test$od_max[data_test$plasmid_type == pair[2]]
        )
        data.frame(
            DHFR_layout1 = unique(data_test$DHFR_layout[data_test$plasmid_type == pair[1]]),
            DHFR_layout2 = unique(data_test$DHFR_layout[data_test$plasmid_type == pair[2]]),
            meaning1 = unique(data_test$meaning[data_test$plasmid_type == pair[1]]),
            meaning2 = unique(data_test$meaning[data_test$plasmid_type == pair[2]]),
            p_value = test$p.value,
            cohen_d = cd$estimate
        )
    })) %>%
        as.data.frame() %>%
        mutate(
            p_value = as.numeric(p_value),
            cohen_d = as.numeric(cohen_d),
            # DHFR_layout1 = factor(DHFR_layout1, level = rev(unique(DHFR_layout1)))
            # DHFR_layout2 = factor(DHFR_layout2, level = rev(unique(DHFR_layout1)))
        )
    # set color bar
    color_pvalue <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")[1:6]))(255)
    color_cohend <- colorRampPalette((brewer.pal(11, "RdBu")[2:10]))(255)
    p_t_odmax <- ggplot(data = t_odmax, aes(x = DHFR_layout2, y = DHFR_layout1)) +
        facet_grid(meaning1 ~ meaning2, scales = "fixed") +
        coord_fixed(ratio = 1) +
        geom_tile(aes(fill = p_value)) +
        geom_text(aes(label = format(p_value, scientific = TRUE, digits = 3))) +
        scale_fill_gradientn(colours = color_pvalue) +
        labs(title = expression("T-test: " * H[0] * ":" * OD600[max[DHFR_layout1]] - OD600[max[DHFR_layout2]] <= 0)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    ggsave(
        plot = p_t_odmax,
        filename = paste("figure_", date, "/t_test_odmax.png", sep = ""),
        width = 11, height = 9
    )
    print("t_test_odmax finish.")
    p_cd_odmax <- ggplot(data = t_odmax, aes(x = DHFR_layout2, y = DHFR_layout1)) +
        facet_grid(meaning1 ~ meaning2, scales = "fixed") +
        coord_fixed(ratio = 1) +
        geom_tile(aes(fill = cohen_d)) +
        geom_text(aes(label = round(cohen_d, 2))) +
        scale_fill_gradientn(colours = color_cohend) +
        labs(title = expression("Cohen's d:" * OD600[max[DHFR_layout1]] - OD600[max[DHFR_layout2]])) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    ggsave(
        plot = p_cd_odmax,
        filename = paste("figure_", date, "/cohen_d_odmax.png", sep = ""),
        width = 11, height = 9
    )
    print("cohen_d_odmax finish.")
    #---- growth_rate----
    t_gr <- do.call(rbind, lapply(all_pairs, function(pair) {
        test <- t.test(
            data_test$growth_rate[data_test$plasmid_type == pair[1]],
            data_test$growth_rate[data_test$plasmid_type == pair[2]],
            alternative = "greater"
        )
        cd <- cohen.d(
            data_test$growth_rate[data_test$plasmid_type == pair[1]],
            data_test$growth_rate[data_test$plasmid_type == pair[2]]
        )
        data.frame(
            DHFR_layout1 = unique(data_test$DHFR_layout[data_test$plasmid_type == pair[1]]),
            DHFR_layout2 = unique(data_test$DHFR_layout[data_test$plasmid_type == pair[2]]),
            meaning1 = unique(data_test$meaning[data_test$plasmid_type == pair[1]]),
            meaning2 = unique(data_test$meaning[data_test$plasmid_type == pair[2]]),
            critical_value = test$conf.int[1],
            p_value = test$p.value,
            cohen_d = cd$estimate
        )
    })) %>%
        as.data.frame() %>%
        mutate(
            p_value = as.numeric(p_value),
            cohen_d = as.numeric(cohen_d),
            # DHFR_layout1 = factor(DHFR_layout1, levels = rev(unique(DHFR_layout1))),
            # DHFR_layout2 = factor(DHFR_layout2, levels = rev(unique(DHFR_layout1)))
        )
    p_t_gr <- ggplot(data = t_gr, aes(x = DHFR_layout2, y = DHFR_layout1)) +
        facet_grid(meaning1 ~ meaning2, scales = "fixed") +
        coord_fixed(ratio = 1) +
        geom_tile(aes(fill = p_value)) +
        geom_text(aes(label = format(p_value, scientific = TRUE, digits = 3))) +
        scale_fill_gradientn(colours = color_pvalue) +
        labs(title = expression("T-test: " * H[0] * ":" * Growth_rate[DHFR_layout1] - Growth_rate[DHFR_layout2] <= 0)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    ggsave(
        plot = p_t_gr,
        filename = paste("figure_", date, "/t_test_gr.png", sep = ""),
        width = 11, height = 9
    )
    print("t_test_gr finish.")
    p_cd_odmax <- ggplot(data = t_gr, aes(x = DHFR_layout2, y = DHFR_layout1)) +
        facet_grid(meaning1 ~ meaning2, scales = "fixed") +
        coord_fixed(ratio = 1) +
        geom_tile(aes(fill = cohen_d)) +
        geom_text(aes(label = round(cohen_d, 2))) +
        scale_fill_gradientn(colours = color_cohend) +
        labs(title = expression("Cohen's d:" * Growth_rate[DHFR_layout1] - Growth_rate[DHFR_layout2])) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
    ggsave(
        plot = p_cd_odmax,
        filename = paste("figure_", date, "/cohen_d_gr.png", sep = ""),
        width = 11, height = 9
    )
    print("cohen_d_gr finish.")
}

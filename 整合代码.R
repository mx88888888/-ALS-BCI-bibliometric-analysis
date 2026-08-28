# ================================================================
# ALS-BCI 文献计量学完整分析脚本（全新重写版）
# 数据来源：WoS核心合集 2011-2026
# 输出：总体概览 + 合作网络 + 侵入/非侵入对比
# 重点修复：作者拆分（strsplit）、关键词合并（强制大写+映射）、机构标准化
# ================================================================

# ================================================================
# 第0部分：环境配置
# ================================================================

rm(list = ls())

# 加载包（如未安装请先 install.packages()）
library(bibliometrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)
library(stringr)

# 设置工作目录（修改为你的实际路径）
setwd("E:/脑机接口对比文献计量学/2.代码/数据收集")

cat("============================================================\n")
cat("     ALS-BCI 文献计量学完整分析（全新版）\n")
cat("============================================================\n\n")

# ================================================================
# 第1部分：数据加载
# ================================================================

cat("【1】加载数据...\n")

file_all <- "all.txt"
file_inv <- "invasive.txt"
file_non <- "non-invasive.txt"

M_all <- convert2df(file_all, dbsource = "wos", format = "plaintext")
M_inv <- convert2df(file_inv, dbsource = "wos", format = "plaintext")
M_non <- convert2df(file_non, dbsource = "wos", format = "plaintext")

cat("  总体文献:", nrow(M_all), "篇\n")
cat("  侵入式:", nrow(M_inv), "篇\n")
cat("  非侵入式:", nrow(M_non), "篇\n\n")

# ================================================================
# 第2部分：关键词同义词合并（强制大写 + 扩展映射）
# ================================================================

cat("【2】关键词同义词合并...\n")

# 同义词映射表（全部目标词转为大写）
synonym_map <- c(
  # BCI系列 → BRAIN-COMPUTER INTERFACE
  "BCI" = "BRAIN-COMPUTER INTERFACE",
  "BMI" = "BRAIN-COMPUTER INTERFACE",
  "BRAIN-COMPUTER INTERFACE (BCI)" = "BRAIN-COMPUTER INTERFACE",
  "BRAIN-COMPUTER INTERFACES" = "BRAIN-COMPUTER INTERFACE",
  "BRAIN-COMPUTER INTERFACES (BCI)" = "BRAIN-COMPUTER INTERFACE",
  "BRAIN-COMPUTER-INTERFACES" = "BRAIN-COMPUTER INTERFACE",
  "BRAIN-COMPUTER-INTERFACE" = "BRAIN-COMPUTER INTERFACE",
  "BRAIN-MACHINE INTERFACE" = "BRAIN-COMPUTER INTERFACE",
  "BRAIN MACHINE INTERFACE" = "BRAIN-COMPUTER INTERFACE",
  "BRAIN COMPUTER INTERFACE" = "BRAIN-COMPUTER INTERFACE",
  "BRAIN-MACHINE INTERFACE (BMI)" = "BRAIN-COMPUTER INTERFACE",
  "BRAINCOMPUTER INTERFACE" = "BRAIN-COMPUTER INTERFACE",
  
  # ALS系列 → AMYOTROPHIC LATERAL SCLEROSIS
  "ALS" = "AMYOTROPHIC LATERAL SCLEROSIS",
  "AMYOTROPHIC LATERAL SCLEROSIS (ALS)" = "AMYOTROPHIC LATERAL SCLEROSIS",
  "AMYOTROPHIC-LATERAL-SCLEROSIS" = "AMYOTROPHIC LATERAL SCLEROSIS",
  "LOU GEHRIG DISEASE" = "AMYOTROPHIC LATERAL SCLEROSIS",
  
  # EEG系列 → ELECTROENCEPHALOGRAPHY
  "EEG" = "ELECTROENCEPHALOGRAPHY",
  "ELECTROENCEPHALOGRAPHY (EEG)" = "ELECTROENCEPHALOGRAPHY",
  "ELECTROENCEPHALOGRAM" = "ELECTROENCEPHALOGRAPHY",
  
  # ECoG系列 → ELECTROCORTICOGRAPHY
  "ECOG" = "ELECTROCORTICOGRAPHY",
  "ELECTROCORTICOGRAPHY (ECOG)" = "ELECTROCORTICOGRAPHY",
  
  # P300系列 → P300
  "P300 SPELLER" = "P300",
  "P3" = "P300",
  
  # SSVEP系列 → STEADY-STATE VISUAL EVOKED POTENTIAL
  "SSVEP" = "STEADY-STATE VISUAL EVOKED POTENTIAL",
  "STEADY-STATE VISUAL EVOKED POTENTIAL (SSVEP)" = "STEADY-STATE VISUAL EVOKED POTENTIAL",
  "STEADY-STATE VISUAL EVOKED POTENTIALS" = "STEADY-STATE VISUAL EVOKED POTENTIAL"
)

# 合并函数
merge_keywords <- function(kw_str) {
  if (is.na(kw_str) || kw_str == "") return("")
  kw_list <- unlist(strsplit(kw_str, "; "))
  kw_list <- trimws(kw_list)
  kw_list <- kw_list[nchar(kw_list) > 0]
  kw_list <- toupper(kw_list)  # 强制大写
  for (i in seq_along(kw_list)) {
    if (kw_list[i] %in% names(synonym_map)) {
      kw_list[i] <- synonym_map[kw_list[i]]
    }
  }
  kw_list <- unique(kw_list)
  return(paste(kw_list, collapse = "; "))
}

# 应用到三个数据集
M_all$DE <- sapply(M_all$DE, merge_keywords)
M_inv$DE <- sapply(M_inv$DE, merge_keywords)
M_non$DE <- sapply(M_non$DE, merge_keywords)

# 验证合并效果
cat("  验证合并效果:\n")
for (df_name in c("M_all", "M_inv", "M_non")) {
  df <- get(df_name)
  kw_test <- unlist(strsplit(df$DE, "; "))
  kw_test <- kw_test[nchar(kw_test) > 0]
  cat("    ", df_name, "关键词总数:", length(kw_test),
      " | BCI残留:", sum(kw_test == "BCI"),
      " | EEG残留:", sum(kw_test == "EEG"),
      " | ALS残留:", sum(kw_test == "ALS"), "\n")
}
cat("  ✅ 关键词合并完成\n\n")

# ================================================================
# 第3部分：机构名称标准化
# ================================================================

cat("【3】机构名称标准化...\n")

# 提取AU_UN字段
if (!("AU_UN" %in% colnames(M_all))) {
  M_all <- metaTagExtraction(M_all, Field = "AU_UN", sep = ";")
}
if (!("AU_UN" %in% colnames(M_inv))) {
  M_inv <- metaTagExtraction(M_inv, Field = "AU_UN", sep = ";")
}
if (!("AU_UN" %in% colnames(M_non))) {
  M_non <- metaTagExtraction(M_non, Field = "AU_UN", sep = ";")
}

# 机构标准化函数
# ---- 机构标准化函数（最终完整版：加州大学全部合并为系统） ----
standardize_affiliation <- function(aff_str) {
  if (is.na(aff_str) || aff_str == "") return("")
  
  # 统一用 ";" 拆分
  inst_list <- unlist(strsplit(aff_str, ";"))
  inst_list <- trimws(inst_list)
  inst_list <- toupper(inst_list)
  inst_list <- inst_list[nchar(inst_list) > 0]
  
  if (length(inst_list) == 0) return("")
  
  for (i in seq_along(inst_list)) {
    x <- inst_list[i]
    
    # === 1. 蒂宾根大学（含 EBERHARD KARLS 变体）===
    if (grepl("TUBINGEN|TUEBINGEN|EBERHARD KARLS", x)) {
      inst_list[i] <- "UNIVERSITY OF TUBINGEN"
      next
    }
    
    # === 2. 宾夕法尼亚州立大学 ===
    if (grepl("PENNSYLVANIA COMMONWEALTH|PENN STATE|PENNSYLVANIA STATE", x)) {
      inst_list[i] <- "PENNSYLVANIA STATE UNIVERSITY"
      next
    }
    
    # === 3. 约翰霍普金斯大学 ===
    if (grepl("JOHNS HOPKINS", x)) {
      inst_list[i] <- "JOHNS HOPKINS UNIVERSITY"
      next
    }
    
    # === 4. 斯坦福大学 ===
    if (grepl("STANFORD", x)) {
      inst_list[i] <- "STANFORD UNIVERSITY"
      next
    }
    
    # === 5. 哈佛大学 ===
    if (grepl("HARVARD", x)) {
      inst_list[i] <- "HARVARD UNIVERSITY"
      next
    }
    
    # === 6. 麻省总医院 ===
    if (grepl("MASSACHUSETTS GENERAL", x)) {
      inst_list[i] <- "MASSACHUSETTS GENERAL HOSPITAL"
      next
    }
    
    # === 7. 布朗大学 ===
    if (grepl("BROWN", x) && !grepl("UNIVERSITY OF CALIFORNIA", x)) {
      inst_list[i] <- "BROWN UNIVERSITY"
      next
    }
    
    # === 8. 墨尔本大学 ===
    if (grepl("MELBOURNE", x)) {
      inst_list[i] <- "UNIVERSITY OF MELBOURNE"
      next
    }
    
    # === 9. 匹兹堡大学 ===
    if (grepl("UNIVERSITY OF PITTSBURGH", x)) {
      inst_list[i] <- "UNIVERSITY OF PITTSBURGH"
      next
    }
    
    # === 10. 密歇根大学 ===
    if (grepl("UNIVERSITY OF MICHIGAN", x)) {
      inst_list[i] <- "UNIVERSITY OF MICHIGAN"
      next
    }
    
    # === 11. 乌得勒支大学 ===
    if (grepl("UTRECHT", x)) {
      inst_list[i] <- "UTRECHT UNIVERSITY"
      next
    }
    
    # === 12. 罗马大学 ===
    if (grepl("SAPIENZA", x)) {
      inst_list[i] <- "SAPIENZA UNIVERSITY ROME"
      next
    }
    
    # === 13. 弗莱堡大学 ===
    if (grepl("FREIBURG", x)) {
      inst_list[i] <- "UNIVERSITY OF FREIBURG"
      next
    }
    
    # === 14. 大阪大学 ===
    if (grepl("OSAKA", x)) {
      inst_list[i] <- "UNIVERSITY OF OSAKA"
      next
    }
    
    # === 15. 牛津大学 ===
    if (grepl("UNIVERSITY OF OXFORD|OXFORD UNIV", x)) {
      inst_list[i] <- "UNIVERSITY OF OXFORD"
      next
    }
    
    # === 16. 科英布拉大学 ===
    if (grepl("COIMBRA", x)) {
      inst_list[i] <- "UNIVERSIDADE DE COIMBRA"
      next
    }
    
    # === 17. 维尔茨堡大学 ===
    if (grepl("WURZBURG", x)) {
      inst_list[i] <- "UNIVERSITY OF WURZBURG"
      next
    }
    
    # === 18. 卡内基梅隆大学 ===
    if (grepl("CARNEGIE MELLON", x)) {
      inst_list[i] <- "CARNEGIE MELLON UNIVERSITY"
      next
    }
    
    # === 19. 马普所 ===
    if (grepl("MAX PLANCK", x)) {
      inst_list[i] <- "MAX PLANCK SOCIETY"
      next
    }
    
    # === 20. 纽约州立大学 ===
    if (grepl("STATE UNIVERSITY OF NEW YORK|SUNY", x)) {
      inst_list[i] <- "STATE UNIVERSITY OF NEW YORK"
      next
    }
    
    # === 21. 奥尔堡大学 ===
    if (grepl("AALBORG", x)) {
      inst_list[i] <- "AALBORG UNIVERSITY"
      next
    }
    
    # === 22. 奥胡斯大学 ===
    if (grepl("AARHUS", x)) {
      inst_list[i] <- "AARHUS UNIVERSITY"
      next
    }
    
    # === 23. 千叶大学 ===
    if (grepl("CHIBA", x)) {
      inst_list[i] <- "CHIBA UNIVERSITY"
      next
    }
    
    # === 24. 帕多瓦大学 ===
    if (grepl("PADUA|PADOVA", x)) {
      inst_list[i] <- "UNIVERSITY OF PADUA"
      next
    }
    
    # === 25. 马拉加大学 ===
    if (grepl("MALAGA", x)) {
      inst_list[i] <- "UNIVERSIDAD DE MALAGA"
      next
    }
    
    # === 26. 蒙特雷理工学院 ===
    if (grepl("TECNOLOGICO DE MONTERREY|MONTERREY", x) && !grepl("TEXAS", x)) {
      inst_list[i] <- "TECNOLOGICO DE MONTERREY"
      next
    }
    
    # === 27. 北航 ===
    if (grepl("BEIHANG", x)) {
      inst_list[i] <- "BEIHANG UNIVERSITY"
      next
    }
    
    # === 28. 北京工业大学 ===
    if (grepl("BEIJING UNIVERSITY OF TECHNOLOGY", x)) {
      inst_list[i] <- "BEIJING UNIVERSITY OF TECHNOLOGY"
      next
    }
    
    # === 29. 首都医科大学 ===
    if (grepl("CAPITAL MEDICAL", x)) {
      inst_list[i] <- "CAPITAL MEDICAL UNIVERSITY"
      next
    }
    
    # === 30. 天津工业大学 ===
    if (grepl("TIANJIN UNIVERSITY OF TECHNOLOGY", x)) {
      inst_list[i] <- "TIANJIN UNIVERSITY OF TECHNOLOGY"
      next
    }
    
    # === 31. 帝国理工 ===
    if (grepl("IMPERIAL COLLEGE", x)) {
      inst_list[i] <- "IMPERIAL COLLEGE LONDON"
      next
    }
    
    # === 32. 伦敦大学学院 ===
    if (grepl("UNIVERSITY COLLEGE LONDON|UCL", x) && !grepl("HOSPITAL", x)) {
      inst_list[i] <- "UNIVERSITY COLLEGE LONDON"
      next
    }
    
    # === 33. 加州大学系统（★ 全部合并为系统名 ★）===
    if (grepl("UNIVERSITY OF CALIFORNIA|UC DAVIS|UCLA|UCSF|UC BERKELEY|UC IRVINE|UC SAN DIEGO|UC SAN FRANCISCO", x)) {
      inst_list[i] <- "UNIVERSITY OF CALIFORNIA SYSTEM"
      next
    }
    
    # === 34. CNRS ===
    if (grepl("CENTRE NATIONAL DE LA RECHERCHE|CNRS", x)) {
      inst_list[i] <- "CNRS"
      next
    }
    
    # === 35. INSERM ===
    if (grepl("INSERM", x)) {
      inst_list[i] <- "INSERM"
      next
    }
    
    # === 36. IRCCS 保留原名 ===
    if (grepl("IRCCS", x)) {
      inst_list[i] <- x
      next
    }
    
    # === 37. VA 系统 ===
    if (grepl("VETERANS AFFAIRS|VHA", x)) {
      inst_list[i] <- "US DEPARTMENT OF VETERANS AFFAIRS"
      next
    }
    
    # === 38. 通用清理 ===
    x <- gsub("^UNIV OF ", "UNIVERSITY OF ", x)
    x <- gsub("^UNIV ", "UNIVERSITY OF ", x)
    x <- gsub(" SYSTEM$", "", x)
    
    inst_list[i] <- x
  }
  
  # 去重
  inst_list <- unique(inst_list)
  
  # 按字母排序
  inst_list <- sort(inst_list)
  
  return(paste(inst_list, collapse = "; "))
}
M_all$AU_UN_Clean <- sapply(M_all$AU_UN, standardize_affiliation)
M_inv$AU_UN_Clean <- sapply(M_inv$AU_UN, standardize_affiliation)
M_non$AU_UN_Clean <- sapply(M_non$AU_UN, standardize_affiliation)
# ---- 从 M_all 中提取全部机构名称 ----
all_aff <- unlist(strsplit(M_all$AU_UN, ";"))
all_aff <- trimws(all_aff)
all_aff <- toupper(all_aff)
all_aff <- all_aff[nchar(all_aff) > 0]

# 统计每个名称出现的次数
aff_freq <- as.data.frame(sort(table(all_aff), decreasing = TRUE))
names(aff_freq) <- c("Raw_Affiliation", "Freq")

# 查看前50条
head(aff_freq, 50)

# 保存为CSV，发给我
write.csv(aff_freq, "All_Raw_Affiliations_Full.csv", row.names = FALSE)
cat("  ✅ 机构标准化完成\n\n")

# ================================================================
# 第4部分：总体概览分析
# ================================================================

cat("【4】总体概览分析\n")

# ---------- 4.1 发文量趋势 ----------
cat("  4.1 发文量趋势...\n")

py_all <- table(M_all$PY)
py_df <- data.frame(Year = as.integer(names(py_all)), 
                    Count = as.integer(py_all))

p1 <- ggplot(py_df, aes(x = Year, y = Count)) +
  geom_line(color = "#2C3E50", size = 1.2) +
  geom_point(color = "#E74C3C", size = 3) +
  labs(title = "ALS-BCI Annual Publication Output (2011-2026)",
       x = "Year", y = "Number of Publications") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Figure1_Overall_PY_Trend.png", p1, width = 10, height = 6, dpi = 300)
write.csv(py_df, "Table_Overall_PY.csv", row.names = FALSE)
cat("    ✅ Figure1_Overall_PY_Trend.png\n")

# ---------- 4.2 核心作者（普赖斯定律）----------
cat("  4.2 核心作者统计...\n")

# 拆分作者
auth_list <- strsplit(M_all$AU, ";")
all_authors <- unlist(auth_list)
all_authors <- trimws(all_authors)
all_authors <- all_authors[all_authors != ""]

author_freq <- as.data.frame(table(all_authors), stringsAsFactors = FALSE)
names(author_freq) <- c("Author", "Documents")
author_freq <- author_freq[order(-author_freq$Documents), ]

# 普赖斯定律
Nmax <- max(author_freq$Documents)
m_threshold <- ceiling(0.749 * sqrt(Nmax))
core_authors <- author_freq[author_freq$Documents >= m_threshold, ]
total_docs <- nrow(M_all)

# ★ 修正：计算“包含至少一位核心作者的文章数” ★
core_names <- core_authors$Author
has_core <- sapply(M_all$AU, function(x) {
  if (is.na(x) || x == "") return(FALSE)
  authors_in_paper <- trimws(unlist(strsplit(x, ";")))
  any(authors_in_paper %in% core_names)
})
papers_with_core <- sum(has_core)
core_ratio <- papers_with_core / total_docs

cat("    总文章数:", total_docs, "\n")
cat("    唯一作者数:", nrow(author_freq), "\n")
cat("    最高产作者发文:", Nmax, "\n")
cat("    核心作者阈值(发文≥):", m_threshold, "\n")
cat("    核心作者人数:", nrow(core_authors), "\n")
cat("    包含核心作者的文章数:", papers_with_core, "/", total_docs, "\n")
cat("    核心作者发文占比:", round(core_ratio * 100, 2), "%\n")

# 前10位核心作者
cat("\n    前10位核心作者:\n")
print(head(core_authors, 10))

write.csv(author_freq, "Table_All_Authors_Freq.csv", row.names = FALSE)
write.csv(core_authors, "Table_Core_Authors.csv", row.names = FALSE)
cat("    ✅ 已保存核心作者表\n")

# ---------- 4.3 前20作者多指标（h指数等）----------
cat("  4.3 前20作者多指标...\n")

# 构建作者-论文数据用于计算h指数
df_temp <- data.frame(
  AU = M_all$AU,
  TC = ifelse(is.na(M_all$TC), 0, M_all$TC),
  stringsAsFactors = FALSE
)

# 拆分作者
auth_list2 <- strsplit(df_temp$AU, ";")
# 展开为长格式
author_rows <- data.frame(
  Author = unlist(auth_list2),
  TC = rep(df_temp$TC, times = sapply(auth_list2, length)),
  stringsAsFactors = FALSE
)
author_rows$Author <- trimws(author_rows$Author)
author_rows <- author_rows[author_rows$Author != "", ]

# 计算指标
author_metrics <- author_rows %>%
  group_by(Author) %>%
  summarise(
    Documents = n(),
    Total_Cites = sum(TC, na.rm = TRUE),
    TC_list = list(TC),
    .groups = "drop"
  ) %>%
  mutate(Avg_Cites = round(Total_Cites / Documents, 2)) %>%
  rowwise() %>%
  mutate(
    h_index = {
      cites <- sort(unlist(TC_list), decreasing = TRUE)
      h <- 0
      for (i in seq_along(cites)) {
        if (cites[i] >= i) h <- i else break
      }
      h
    },
    g_index = {
      cites <- sort(unlist(TC_list), decreasing = TRUE)
      total <- 0
      g <- 0
      for (i in seq_along(cites)) {
        total <- total + cites[i]
        if (total >= i^2) g <- i else break
      }
      g
    },
    i10_index = sum(unlist(TC_list) >= 10, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(desc(Documents), desc(Total_Cites)) %>%
  mutate(Rank = row_number()) %>%
  select(Rank, Author, Documents, Total_Cites, Avg_Cites, h_index, g_index, i10_index)

write.csv(head(author_metrics, 20), "Table_Top20_Authors_Metrics.csv", row.names = FALSE)
cat("    ✅ 已保存前20作者指标\n")

# ---------- 4.4 期刊分布 ----------
cat("  4.4 期刊分布...\n")

journal_freq <- as.data.frame(sort(table(M_all$SO), decreasing = TRUE))
names(journal_freq) <- c("Journal", "N")
journal_freq$Percentage <- round(journal_freq$N / nrow(M_all) * 100, 2)
write.csv(head(journal_freq, 20), "Table_Journals_Top20.csv", row.names = FALSE)
cat("    ✅ 已保存期刊分布\n")

# ---------- 4.5 机构分布（使用清洗后字段）----------
cat("  4.5 机构分布...\n")

inst_all <- unlist(strsplit(M_all$AU_UN_Clean, "; "))
inst_all <- inst_all[nchar(inst_all) > 0]
inst_freq <- as.data.frame(sort(table(inst_all), decreasing = TRUE))
names(inst_freq) <- c("Affiliation", "N")
inst_freq$Percentage <- round(inst_freq$N / nrow(M_all) * 100, 2)
write.csv(head(inst_freq, 20), "Table_Institutions_Top20.csv", row.names = FALSE)
cat("    ✅ 已保存机构分布\n")


# ---------- 4.6 高被引文献 ----------
cat("  4.6 高被引文献...\n")

top_papers <- M_all[order(-M_all$TC), c("TI", "AU", "PY", "SO", "TC")]
names(top_papers) <- c("Title", "Authors", "Year", "Journal", "Total_Citations")
write.csv(head(top_papers, 15), "Table_Top_Cited_Papers.csv", row.names = FALSE)
cat("    ✅ 已保存高被引文献\n")

# ---------- 4.7 总体关键词高频统计 ----------
cat("  4.7 总体关键词统计...\n")

kw_all <- unlist(strsplit(M_all$DE, "; "))
kw_all <- kw_all[nchar(kw_all) > 0]
kw_freq <- as.data.frame(sort(table(kw_all), decreasing = TRUE))
names(kw_freq) <- c("Keyword", "Freq")
write.csv(head(kw_freq, 20), "Table_Keywords_Overall_Top20.csv", row.names = FALSE)
cat("    ✅ 已保存总体关键词\n\n")

# ================================================================
# 第5部分：侵入式 vs 非侵入式对比分析
# ================================================================

cat("【5】侵入式 vs 非侵入式对比分析\n")

# ---------- 5.1 发文量与累计趋势（含拐点量化）----------
# ---------- 5.1 发文量与累计趋势（含拐点量化）----------
cat("  5.1 发文量与拐点分析...\n")

py_inv <- table(M_inv$PY)
py_non <- table(M_non$PY)
all_years <- sort(union(names(py_inv), names(py_non)))

py_comp <- data.frame(
  Year = as.integer(all_years),
  Invasive = as.integer(py_inv[all_years]),
  Noninvasive = as.integer(py_non[all_years])
)
py_comp[is.na(py_comp)] <- 0
py_comp$Total <- py_comp$Invasive + py_comp$Noninvasive
py_comp$Cum_Invasive <- cumsum(py_comp$Invasive)
py_comp$Cum_Noninvasive <- cumsum(py_comp$Noninvasive)

# ==== 拐点斜率计算（侵入式 + 非侵入式）====

# 侵入式
inv_pre <- py_comp[py_comp$Year >= 2011 & py_comp$Year <= 2020, ]
inv_post <- py_comp[py_comp$Year >= 2021& py_comp$Year <= 2024, ]
slope_inv_pre <- if(nrow(inv_pre) > 1) coef(lm(Invasive ~ Year, data = inv_pre))[2] else NA
slope_inv_post <- if(nrow(inv_post) > 1) coef(lm(Invasive ~ Year, data = inv_post))[2] else NA

# 非侵入式
non_pre <- py_comp[py_comp$Year >= 2011 & py_comp$Year <= 2020, ]
non_post <- py_comp[py_comp$Year >= 2021& py_comp$Year <= 2024, ]
slope_non_pre <- if(nrow(non_pre) > 1) coef(lm(Noninvasive ~ Year, data = non_pre))[2] else NA
slope_non_post <- if(nrow(non_post) > 1) coef(lm(Noninvasive ~ Year, data = non_post))[2] else NA

# 输出
cat("    ===== 侵入式 =====\n")
cat("      2011-2020 斜率:", round(slope_inv_pre, 3), "\n")
cat("      2021-2024 斜率:", round(slope_inv_post, 3), "\n")
cat("      增长倍数:", round(slope_inv_post / slope_inv_pre, 2), "x\n")
cat("    ===== 非侵入式 =====\n")
cat("      2011-2020 斜率:", round(slope_non_pre, 3), "\n")
cat("      2021-2024 斜率:", round(slope_non_post, 3), "\n")
cat("      增长倍数:", round(slope_non_post / slope_non_pre, 2), "x\n")

# 保存拐点数据
slope_df <- data.frame(
  Type = c("Invasive", "Invasive", "Non-invasive", "Non-invasive"),
  Period = c("2011-2020", "2021-2024", "2011-2020", "2021-2024"),
  Slope = c(slope_inv_pre, slope_inv_post, slope_non_pre, slope_non_post)
)
write.csv(slope_df, "Table_Inflection_Slopes.csv", row.names = FALSE)

write.csv(py_comp, "Table_PY_Comparison.csv", row.names = FALSE)

# 发文量趋势图
py_long <- pivot_longer(py_comp, cols = c(Invasive, Noninvasive), 
                        names_to = "Type", values_to = "Count")

p2 <- ggplot(py_long, aes(x = Year, y = Count, color = Type, group = Type)) +
  geom_line(size = 1.2) + geom_point(size = 2.5) +
  labs(title = "Publication Volume: Invasive vs Non-invasive",
       x = "Year", y = "Number of Publications") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("Invasive" = "#E41A1C", "Noninvasive" = "#377EB8"))
ggsave("Figure2_PY_Comparison.png", p2, width = 10, height = 6, dpi = 300)

# 累计趋势图
py_cum_long <- py_comp %>%
  select(Year, Cum_Invasive, Cum_Noninvasive) %>%
  pivot_longer(cols = c(Cum_Invasive, Cum_Noninvasive), 
               names_to = "Type", values_to = "Cumulative")

p3 <- ggplot(py_cum_long, aes(x = Year, y = Cumulative, color = Type, group = Type)) +
  geom_line(size = 1.2) + geom_point(size = 2) +
  labs(title = "Cumulative Publication Volume", 
       x = "Year", y = "Accumulated Publications") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_manual(values = c("Cum_Invasive" = "#E41A1C", "Cum_Noninvasive" = "#377EB8"),
                     labels = c("Invasive", "Non-invasive"))
ggsave("Figure3_Cumulative_Comparison.png", p3, width = 10, height = 6, dpi = 300)
cat("    ✅ 已保存趋势图和拐点斜率表\n")

# ---------- 5.2 期刊分布对比 ----------
cat("  5.2 期刊分布对比...\n")

j_inv <- sort(table(M_inv$SO), decreasing = TRUE)
j_non <- sort(table(M_non$SO), decreasing = TRUE)
j_df_inv <- data.frame(Journal = names(j_inv), N = as.integer(j_inv))
j_df_non <- data.frame(Journal = names(j_non), N = as.integer(j_non))

write.csv(head(j_df_inv, 15), "Table_Journals_Invasive.csv", row.names = FALSE)
write.csv(head(j_df_non, 15), "Table_Journals_Noninvasive.csv", row.names = FALSE)

# 期刊对比图
j_top_inv <- head(j_df_inv, 10); j_top_inv$Group <- "Invasive"
j_top_non <- head(j_df_non, 10); j_top_non$Group <- "Non-invasive"
j_comb <- rbind(j_top_inv, j_top_non)

p4 <- ggplot(j_comb, aes(x = reorder(Journal, N), y = N, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Top 10 Journals: Invasive vs Non-invasive", 
       x = "", y = "Publications") +
  theme_minimal() + theme(legend.position = "bottom") +
  scale_fill_manual(values = c("Invasive" = "#E41A1C", "Non-invasive" = "#377EB8"))
ggsave("Figure4_Journal_Comparison.png", p4, width = 10, height = 7, dpi = 300)
cat("    ✅ 已保存期刊对比图\n")

# ---------- 5.3 国家分布对比 ----------
cat("  5.3 国家分布对比...\n")

# 提取国家字段
if (!("AU_CO" %in% colnames(M_inv))) {
  M_inv <- metaTagExtraction(M_inv, "AU_CO", sep = ";")
}
if (!("AU_CO" %in% colnames(M_non))) {
  M_non <- metaTagExtraction(M_non, "AU_CO", sep = ";")
}

c_inv_raw <- strsplit(M_inv$AU_CO, ";")
c_inv <- unlist(c_inv_raw); c_inv <- trimws(c_inv); c_inv <- c_inv[nchar(c_inv) > 0]
c_non_raw <- strsplit(M_non$AU_CO, ";")
c_non <- unlist(c_non_raw); c_non <- trimws(c_non); c_non <- c_non[nchar(c_non) > 0]

c_df_inv <- as.data.frame(sort(table(c_inv), decreasing = TRUE))
names(c_df_inv) <- c("Country", "Freq")
c_df_non <- as.data.frame(sort(table(c_non), decreasing = TRUE))
names(c_df_non) <- c("Country", "Freq")

write.csv(head(c_df_inv, 15), "Table_Countries_Invasive.csv", row.names = FALSE)
write.csv(head(c_df_non, 15), "Table_Countries_Noninvasive.csv", row.names = FALSE)
cat("    ✅ 已保存国家分布\n")

# ---------- 5.4 作者合作网络密度 ----------
cat("  5.4 作者合作网络...\n")

auth_inv <- biblioNetwork(M_inv, analysis = "collaboration", 
                          network = "authors", sep = ";")
auth_non <- biblioNetwork(M_non, analysis = "collaboration", 
                          network = "authors", sep = ";")

g_auth_inv <- graph_from_adjacency_matrix(auth_inv, weighted = TRUE, mode = "undirected")
g_auth_non <- graph_from_adjacency_matrix(auth_non, weighted = TRUE, mode = "undirected")

auth_net_stats <- data.frame(
  Type = c("Invasive", "Non-invasive"),
  Nodes = c(vcount(g_auth_inv), vcount(g_auth_non)),
  Edges = c(ecount(g_auth_inv), ecount(g_auth_non)),
  Density = c(edge_density(g_auth_inv), edge_density(g_auth_non))
)
write.csv(auth_net_stats, "Table_Auth_Network_Stats.csv", row.names = FALSE)

# 绘制网络图
png("Figure5_Auth_Network_Invasive.png", width = 1200, height = 1000, res = 150)
networkPlot(auth_inv, n = 15, Title = "Invasive: Author Collaboration", 
            type = "fruchterman", labelsize = 0.8, size = TRUE)
dev.off()

png("Figure6_Auth_Network_Noninvasive.png", width = 1200, height = 1000, res = 150)
networkPlot(auth_non, n = 15, Title = "Non-invasive: Author Collaboration", 
            type = "fruchterman", labelsize = 0.8, size = TRUE)
dev.off()
cat("    ✅ 已保存作者合作网络\n")

# ---------- 5.5 期刊共被引网络密度 ----------
cat("  5.5 期刊共被引网络...\n")

if (!("CR_SO" %in% colnames(M_inv))) {
  M_inv <- metaTagExtraction(M_inv, "CR_SO", sep = ";")
}
if (!("CR_SO" %in% colnames(M_non))) {
  M_non <- metaTagExtraction(M_non, "CR_SO", sep = ";")
}

jcr_inv <- biblioNetwork(M_inv, analysis = "co-citation", 
                         network = "sources", sep = ";")
jcr_non <- biblioNetwork(M_non, analysis = "co-citation", 
                         network = "sources", sep = ";")

g_jcr_inv <- graph_from_adjacency_matrix(jcr_inv, weighted = TRUE, mode = "undirected")
g_jcr_non <- graph_from_adjacency_matrix(jcr_non, weighted = TRUE, mode = "undirected")

jcr_stats <- data.frame(
  Type = c("Invasive", "Non-invasive"),
  Nodes = c(vcount(g_jcr_inv), vcount(g_jcr_non)),
  Edges = c(ecount(g_jcr_inv), ecount(g_jcr_non)),
  Density = c(edge_density(g_jcr_inv), edge_density(g_jcr_non))
)
write.csv(jcr_stats, "Table_JCR_Network_Stats.csv", row.names = FALSE)

png("Figure7_JCR_Network_Invasive.png", width = 1000, height = 800, res = 150)
networkPlot(jcr_inv, n = 20, Title = "Invasive: Journal Co-citation", 
            type = "fruchterman", labelsize = 0.8, size = TRUE)
dev.off()

png("Figure8_JCR_Network_Noninvasive.png", width = 1000, height = 800, res = 150)
networkPlot(jcr_non, n = 20, Title = "Non-invasive: Journal Co-citation", 
            type = "fruchterman", labelsize = 0.8, size = TRUE)
dev.off()
cat("    ✅ 已保存期刊共被引网络\n")

# ---------- 5.6 关键词共现网络 ----------
cat("  5.6 关键词共现网络...\n")

Net_inv <- biblioNetwork(M_inv, analysis = "co-occurrences", 
                         network = "author_keywords", sep = ";")
Net_non <- biblioNetwork(M_non, analysis = "co-occurrences", 
                         network = "author_keywords", sep = ";")

g_kw_inv <- graph_from_adjacency_matrix(Net_inv, weighted = TRUE, mode = "undirected")
g_kw_non <- graph_from_adjacency_matrix(Net_non, weighted = TRUE, mode = "undirected")

kw_net_stats <- data.frame(
  Type = c("Invasive", "Non-invasive"),
  Nodes = c(vcount(g_kw_inv), vcount(g_kw_non)),
  Edges = c(ecount(g_kw_inv), ecount(g_kw_non)),
  Density = c(edge_density(g_kw_inv), edge_density(g_kw_non))
)
write.csv(kw_net_stats, "Table_Keyword_Network_Stats.csv", row.names = FALSE)

png("Figure9_Keyword_Network_Invasive.png", width = 1200, height = 1000, res = 150)
networkPlot(Net_inv, n = 15, Title = "Invasive: Keyword Co-occurrence", 
            type = "fruchterman", labelsize = 0.8, size = TRUE)
dev.off()

png("Figure10_Keyword_Network_Noninvasive.png", width = 1200, height = 1000, res = 150)
networkPlot(Net_non, n = 15, Title = "Non-invasive: Keyword Co-occurrence", 
            type = "fruchterman", labelsize = 0.8, size = TRUE)
dev.off()

# 关键词频次对比
kw_inv <- unlist(strsplit(M_inv$DE, "; ")); kw_inv <- kw_inv[nchar(kw_inv) > 0]
kw_non <- unlist(strsplit(M_non$DE, "; ")); kw_non <- kw_non[nchar(kw_non) > 0]
kw_inv_freq <- as.data.frame(sort(table(kw_inv), decreasing = TRUE))
kw_non_freq <- as.data.frame(sort(table(kw_non), decreasing = TRUE))
names(kw_inv_freq) <- c("Keyword", "Freq")
names(kw_non_freq) <- c("Keyword", "Freq")
write.csv(head(kw_inv_freq, 15), "Table_Keywords_Invasive.csv", row.names = FALSE)
write.csv(head(kw_non_freq, 15), "Table_Keywords_Noninvasive.csv", row.names = FALSE)
cat("    ✅ 已保存关键词共现网络\n")

# ---------- 5.7 EEG-AI关联分析 ----------
cat("  5.7 EEG-AI关联分析...\n")

ai_keywords <- c("DEEP LEARNING", "CONVOLUTIONAL NEURAL", "NEURAL NETWORK", 
                 "MACHINE LEARNING", "ARTIFICIAL INTELLIGENCE", "CNN", "RNN")

eeg_ai_papers <- c()
for (i in 1:nrow(M_all)) {
  kw <- M_all$DE[i]
  if (grepl("ELECTROENCEPHALOGRAPHY", kw) && 
      any(sapply(ai_keywords, function(x) grepl(x, kw, fixed = TRUE)))) {
    eeg_ai_papers <- c(eeg_ai_papers, M_all$PY[i])
  }
}

eeg_ai_df <- as.data.frame(table(eeg_ai_papers))
names(eeg_ai_df) <- c("Year", "Count")
eeg_ai_df$Year <- as.integer(as.character(eeg_ai_df$Year))

cat("    EEG-AI共现文献总数:", length(eeg_ai_papers), "\n")
cat("    2021年前数量:", sum(eeg_ai_df$Count[eeg_ai_df$Year < 2021]), "\n")
cat("    2021年后数量:", sum(eeg_ai_df$Count[eeg_ai_df$Year >= 2021]), "\n")
write.csv(eeg_ai_df, "Table_EEG_AI_Association.csv", row.names = FALSE)
cat("    ✅ 已保存EEG-AI关联数据\n")

# ---------- 5.8 伦理词频检索 ----------
cat("  5.8 伦理/患者偏好词频...\n")

ethics_keywords <- c("ETHIC", "ETHICAL", "INFORMED CONSENT", 
                     "PATIENT PREFERENCE", "AUTONOMY", "BIOETHICS", "ETHICS")
kw_all_upper <- toupper(kw_all)
ethics_hits <- sum(sapply(ethics_keywords, function(kw) {
  sum(grepl(kw, kw_all_upper, fixed = TRUE))
}))

cat("    伦理相关关键词出现总频次:", ethics_hits, "\n")
write.csv(data.frame(Keyword = ethics_keywords, 
                     Frequency = ethics_hits), 
          "Table_Ethics_Frequency.csv", row.names = FALSE)
cat("    ✅ 已保存伦理词频\n\n")

# ================================================================
# 第6部分：主题地图（Thematic Map）
# ================================================================

cat("【6】主题地图（Thematic Map）...\n")
cat("    ⚠️ 此步骤可能耗时较长，请耐心等待...\n")

# 注意：thematicMap 可能因数据量或版本问题报错
# 用 tryCatch 包裹，即使失败也不影响其他输出

tryCatch({
  # 侵入式
  cat("    正在生成侵入式主题地图...\n")
  thematic_inv <- thematicMap(M_inv, field = "DE", n = 15, minfreq = 2, 
                              cluster = "louvain")
  png("Figure11_Thematic_Map_Invasive.png", width = 1000, height = 800, res = 150)
  plot(thematic_inv$map)
  dev.off()
  write.csv(thematic_inv$clusters, "Table_Thematic_Clusters_Invasive.csv", row.names = FALSE)
  cat("    ✅ 侵入式主题地图已保存\n")
}, error = function(e) {
  cat("    ⚠️ 侵入式主题地图失败:", e$message, "\n")
})

tryCatch({
  # 非侵入式
  cat("    正在生成非侵入式主题地图...\n")
  thematic_non <- thematicMap(M_non, field = "DE", n = 15, minfreq = 2,
                              cluster = "louvain")
  png("Figure12_Thematic_Map_Noninvasive.png", width = 1000, height = 800, res = 150)
  plot(thematic_non$map)
  dev.off()
  write.csv(thematic_non$clusters, "Table_Thematic_Clusters_Noninvasive.csv", row.names = FALSE)
  cat("    ✅ 非侵入式主题地图已保存\n")
}, error = function(e) {
  cat("    ⚠️ 非侵入式主题地图失败:", e$message, "\n")
})

#高频词共现矩阵热图
# ================================================================
# 关键词共现热图（侵入式 vs 非侵入式）
# 用于替代 Thematic Map，展示高频关键词之间的共现强度
# ================================================================

library(corrplot)

# ---- 定义一个绘制热图的函数 ----
plot_cooccurrence_heatmap <- function(M_data, title_prefix, n_top = 15, exclude_top2 = TRUE) {
  
  # 1. 提取并合并关键词
  all_kw <- unlist(strsplit(M_data$DE, "; "))
  all_kw <- all_kw[nchar(all_kw) > 0]
  
  # 2. 统计频次
  kw_freq <- sort(table(all_kw), decreasing = TRUE)
  
  # 3. 如果 exclude_top2 = TRUE，移除前两个最高频词（通常是 BCI 和 ALS）
  if (exclude_top2 && length(kw_freq) >= 3) {
    kw_freq <- kw_freq[-(1:2)]  # 移除前两个
  }
  
  # 4. 取前 n_top 个高频词
  top_kw <- head(kw_freq, n_top)
  kw_list <- names(top_kw)
  
  cat("  使用的关键词:", paste(kw_list, collapse = ", "), "\n")
  
  # 5. 构建共现矩阵
  co_occur <- matrix(0, nrow = length(kw_list), ncol = length(kw_list))
  rownames(co_occur) <- colnames(co_occur) <- kw_list
  
  for (i in 1:nrow(M_data)) {
    kw_in_doc <- unlist(strsplit(M_data$DE[i], "; "))
    kw_in_doc <- kw_in_doc[nchar(kw_in_doc) > 0]
    
    # 只保留在 kw_list 中的关键词
    kw_in_doc <- kw_in_doc[kw_in_doc %in% kw_list]
    
    if (length(kw_in_doc) >= 2) {
      for (j in 1:(length(kw_in_doc) - 1)) {
        for (k in (j + 1):length(kw_in_doc)) {
          co_occur[kw_in_doc[j], kw_in_doc[k]] <- co_occur[kw_in_doc[j], kw_in_doc[k]] + 1
          co_occur[kw_in_doc[k], kw_in_doc[j]] <- co_occur[kw_in_doc[k], kw_in_doc[j]] + 1  # 对称填充
        }
      }
    }
  }
  
  # 6. 归一化（按行归一化为比例，突出相对共现强度）
  # 避免除以 0
  row_sums <- rowSums(co_occur)
  row_sums[row_sums == 0] <- 1
  co_occur_norm <- co_occur / row_sums
  
  # 7. 绘制热图
  corrplot(co_occur_norm, 
           method = "shade",           # 使用渐变阴影，比 "color" 更清晰
           type = "full",              # 显示完整矩阵
           tl.cex = 0.9,               # 标签大小
           tl.col = "black",           # 标签颜色
           tl.srt = 90,                # 标签旋转角度（便于长词显示）
           cl.cex = 0.8,               # 图例大小
           is.corr = FALSE,            # 不是相关系数矩阵，不强制归一化到 [-1,1]
           col = colorRampPalette(c("white", "#FFE4B5", "#E41A1C"))(100),  # 自定义颜色
           mar = c(2, 2, 2, 2),        # 边距
           main = paste(title_prefix, "Keyword Co-occurrence Heatmap"))
  
  # 8. 返回共现矩阵供后续使用
  return(list(matrix = co_occur, normalized = co_occur_norm, keywords = kw_list))
}

# ---- 绘制侵入式热图 ----
cat("【绘制侵入式关键词共现热图】\n")
result_inv <- plot_cooccurrence_heatmap(
  M_data = M_inv,
  title_prefix = "Invasive",
  n_top = 15,
  exclude_top2 = TRUE
)
# 保存为高分辨率 PNG
png("Figure_Heatmap_Invasive.png", width = 1200, height = 1200, res = 150)
plot_cooccurrence_heatmap(M_inv, "Invasive", 15, TRUE)
dev.off()
cat("  ✅ 已保存: Figure_Heatmap_Invasive.png\n")

# ---- 绘制非侵入式热图 ----
cat("\n【绘制非侵入式关键词共现热图】\n")
png("Figure_Heatmap_Noninvasive.png", width = 1200, height = 1200, res = 150)
plot_cooccurrence_heatmap(M_non, "Non-invasive", 15, TRUE)
dev.off()
cat("  ✅ 已保存: Figure_Heatmap_Noninvasive.png\n")

# ---- 可选：保存共现矩阵为 CSV ----
write.csv(result_inv$matrix, "Table_Cooccurrence_Matrix_Invasive.csv", row.names = TRUE)
write.csv(result_inv$normalized, "Table_Cooccurrence_Matrix_Invasive_Normalized.csv", row.names = TRUE)

# ================================================================
# 第6部分：关键词共现聚类网络（干净版，无标签 + 表格）
# 说明：网络图仅展示聚类结构（颜色块），不显示文字
#       具体关键词清单由 Table_Cluster_*.csv 提供
# ================================================================

cat("【6】关键词共现聚类网络（干净版）...\n")

# ---- 辅助函数：计算每个节点的频次 ----
get_keyword_freq <- function(M_data) {
  all_kw <- unlist(strsplit(M_data$DE, "; "))
  all_kw <- all_kw[nchar(all_kw) > 0]
  return(table(all_kw))
}

# ---- 6.1 侵入式：干净聚类网络 ----
cat("  绘制侵入式聚类网络（无标签）...\n")

if (exists("Net_inv")) {
  g_inv <- graph_from_adjacency_matrix(Net_inv, weighted = TRUE, mode = "undirected")
  g_inv <- delete.vertices(g_inv, degree(g_inv) == 0)
  
  if (vcount(g_inv) > 2) {
    cl_inv <- cluster_louvain(g_inv)
    
    # 配色方案
    colrs <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
    vertex_colors <- colrs[membership(cl_inv) %% length(colrs) + 1]
    
    # 节点大小统一（不区分频次），让聚类结构更清晰
    png("Figure11_Keyword_Cluster_Invasive_Clean.png", width = 1600, height = 1600, res = 200)
    
    par(mar = c(1, 1, 4, 1))
    
    plot(cl_inv, g_inv,
         layout = layout_with_fr(g_inv, niter = 500),
         vertex.label = "",                    # 【关键】不显示任何文字标签
         vertex.size = 8,                      # 统一节点大小
         vertex.color = vertex_colors,
         vertex.frame.color = "white",
         vertex.frame.width = 0.5,
         edge.width = 0.3,
         edge.color = "gray70",
         main = "Invasive: Keyword Co-occurrence Clusters",
         cex.main = 1.4)
    
    # 添加图例（聚类编号 → 颜色）
    legend("bottomright", 
           legend = 1:length(unique(membership(cl_inv))),
           col = colrs[1:length(unique(membership(cl_inv)))],
           pch = 19, pt.cex = 2,
           title = "Cluster",
           bty = "n",
           cex = 1.2)
    
    dev.off()
    
    # ---- 保存聚类表格（供论文使用） ----
    kw_freq_inv <- get_keyword_freq(M_inv)
    cluster_df_inv <- data.frame(
      Cluster = as.integer(membership(cl_inv)),
      Keyword = names(membership(cl_inv)),
      Freq = as.integer(kw_freq_inv[names(membership(cl_inv))])
    )
    cluster_df_inv$Freq[is.na(cluster_df_inv$Freq)] <- 0
    cluster_df_inv <- cluster_df_inv[order(cluster_df_inv$Cluster, -cluster_df_inv$Freq), ]
    write.csv(cluster_df_inv, "Table_Cluster_Invasive.csv", row.names = FALSE)
    
    cat("    ✅ 已保存: Figure11_Keyword_Cluster_Invasive_Clean.png\n")
    cat("    聚类数:", length(unique(membership(cl_inv))), "\n")
    cat("    节点总数:", vcount(g_inv), "\n")
    cat("    表格已保存: Table_Cluster_Invasive.csv\n")
  } else {
    cat("    ⚠️ 侵入式网络节点过少\n")
  }
} else {
  cat("    ⚠️ 未找到 Net_inv\n")
}

# ---- 6.2 非侵入式：干净聚类网络 ----
cat("\n  绘制非侵入式聚类网络（无标签）...\n")

if (exists("Net_non")) {
  g_non <- graph_from_adjacency_matrix(Net_non, weighted = TRUE, mode = "undirected")
  g_non <- delete.vertices(g_non, degree(g_non) == 0)
  
  if (vcount(g_non) > 2) {
    cl_non <- cluster_louvain(g_non)
    
    colrs <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
    vertex_colors <- colrs[membership(cl_non) %% length(colrs) + 1]
    
    png("Figure12_Keyword_Cluster_Noninvasive_Clean.png", width = 1600, height = 1600, res = 200)
    
    par(mar = c(1, 1, 4, 1))
    
    plot(cl_non, g_non,
         layout = layout_with_fr(g_non, niter = 500),
         vertex.label = "",                    # 【关键】不显示任何文字标签
         vertex.size = 8,
         vertex.color = vertex_colors,
         vertex.frame.color = "white",
         vertex.frame.width = 0.5,
         edge.width = 0.3,
         edge.color = "gray70",
         main = "Non-invasive: Keyword Co-occurrence Clusters",
         cex.main = 1.4)
    
    legend("bottomright", 
           legend = 1:length(unique(membership(cl_non))),
           col = colrs[1:length(unique(membership(cl_non)))],
           pch = 19, pt.cex = 2,
           title = "Cluster",
           bty = "n",
           cex = 1.2)
    
    dev.off()
    
    kw_freq_non <- get_keyword_freq(M_non)
    cluster_df_non <- data.frame(
      Cluster = as.integer(membership(cl_non)),
      Keyword = names(membership(cl_non)),
      Freq = as.integer(kw_freq_non[names(membership(cl_non))])
    )
    cluster_df_non$Freq[is.na(cluster_df_non$Freq)] <- 0
    cluster_df_non <- cluster_df_non[order(cluster_df_non$Cluster, -cluster_df_non$Freq), ]
    write.csv(cluster_df_non, "Table_Cluster_Noninvasive.csv", row.names = FALSE)
    
    cat("    ✅ 已保存: Figure12_Keyword_Cluster_Noninvasive_Clean.png\n")
    cat("    聚类数:", length(unique(membership(cl_non))), "\n")
    cat("    节点总数:", vcount(g_non), "\n")
    cat("    表格已保存: Table_Cluster_Noninvasive.csv\n")
  } else {
    cat("    ⚠️ 非侵入式网络节点过少\n")
  }
} else {
  cat("    ⚠️ 未找到 Net_non\n")
}

cat("\n    ✅ 干净聚类网络 + 表格 全部完成！\n")
# ================================================================
# 第7部分：完成
# ================================================================

cat("\n============================================================\n")
cat("🎉 所有分析已完成！\n\n")
cat("输出文件清单（共约25个文件）：\n")
cat("  图表: Figure1 ~ Figure12\n")
cat("  表格: Table_*.csv\n")
cat("  网络统计: Table_*_Network_Stats.csv\n")
cat("============================================================\n")


# 1. 总体基本信息
cat("=== 总体数据 ===\n")
cat("总文献数:", nrow(M_all), "\n")
cat("总作者数:", length(unique(unlist(strsplit(M_all$AU, ";")))), "\n")
cat("总期刊数:", length(unique(M_all$SO)), "\n")
cat("总机构数:", length(unique(unlist(strsplit(M_all$AU_UN_Clean, ";")))), "\n")
cat("总被引次数:", sum(M_all$TC, na.rm = TRUE), "\n")
cat("时间跨度:", min(M_all$PY), "-", max(M_all$PY), "\n\n")

# 2. 核心作者数据
core_df <- read.csv("Table_Core_Authors.csv")
cat("=== 核心作者 ===\n")
cat("最高产作者发文:", max(core_df$Documents), "\n")
cat("核心作者人数:", nrow(core_df), "\n")
cat("核心作者发文占比:", round(sum(core_df$Documents)/nrow(M_all)*100, 2), "%\n")
cat("前5位核心作者:\n")
print(head(core_df$Author, 5))

# 3. 期刊数据
jour_df <- read.csv("Table_Journals_Top20.csv")
cat("\n=== 前3位期刊 ===\n")
print(head(jour_df, 3))

# 4. 机构数据
inst_df <- read.csv("Table_Institutions_Top20.csv")
cat("\n=== 前3位机构 ===\n")
print(head(inst_df, 3))

# 5. 对比数据
py_comp <- read.csv("Table_PY_Comparison.csv")
cat("\n=== 侵入/非侵入发文 ===\n")
cat("侵入式总发文:", sum(py_comp$Invasive), "\n")
cat("非侵入式总发文:", sum(py_comp$Noninvasive), "\n")

# 6. 被引对比（如果你有计算）
cat("\n=== 被引对比（如果已计算）===\n")
tc_inv <- sum(M_inv$TC, na.rm = TRUE)
tc_non <- sum(M_non$TC, na.rm = TRUE)
cat("侵入式总被引:", tc_inv, " 篇均:", round(tc_inv/nrow(M_inv), 2), "\n")
cat("非侵入式总被引:", tc_non, " 篇均:", round(tc_non/nrow(M_non), 2), "\n")

read.csv("Table_JCR_Network_Stats.csv")
read.csv("Table_Inflection_Slopes.csv")
read.csv("Table_Keyword_Network_Stats.csv")


#机构发文补充
# 验证总机构数
length(unique(unlist(strsplit(M_all$AU_UN_Clean, ";"))))  # 应为394
# 验证发文1篇的机构数
inst_all <- unlist(strsplit(M_all$AU_UN_Clean, ";"))
inst_all <- inst_all[nchar(inst_all) > 0]
inst_freq <- table(inst_all)
sum(inst_freq == 1)  # 应为319

# 补充3.1.2
# ---- 总体期刊共被引网络：聚类 + 图 + 表（修正版） ----
if (!("CR_SO" %in% colnames(M_all))) {
  M_all <- metaTagExtraction(M_all, "CR_SO", sep = ";")
}

# 构建共被引网络
jcr_all <- biblioNetwork(M_all, analysis = "co-citation", 
                         network = "sources", sep = ";")

# 转为igraph
library(igraph)
g_jcr <- graph_from_adjacency_matrix(jcr_all, weighted = TRUE, mode = "undirected")

# ★ 关键：剔除孤立节点（degree=0）
g_jcr <- delete.vertices(g_jcr, degree(g_jcr) == 0)

# Louvain聚类
cl <- cluster_louvain(g_jcr)

# ===== 图：聚类网络图（修正） =====
# 计算节点大小（基于频次）
journal_freq <- sort(table(M_all$SO), decreasing = TRUE)
V(g_jcr)$size <- log(journal_freq[V(g_jcr)$name] + 1) * 3

# ★ 手动指定布局，避免 xlim 报错
set.seed(123)
layout_coords <- layout_with_fr(g_jcr, niter = 500)

png("Figure_JCR_Overall.png", width = 1400, height = 1100, res = 180)
networkPlot(jcr_all, n = 30, 
            Title = "Journal Co-citation Network (Overall)", 
            type = "fruchterman", 
            labelsize = 0.6,
            size = TRUE)
dev.off()

cat("✅ 已生成: Figure_JCR_Overall.png\n")

# ===== 表：聚类成员表 =====
cluster_df <- data.frame(
  Cluster = membership(cl),
  Journal = V(g_jcr)$name
)
write.csv(cluster_df, "Table_JCR_Overall_Clusters.csv", row.names = FALSE)
cat("✅ 已生成: Table_JCR_Overall_Clusters.csv\n")

# 预览聚类结构
cat("\n各聚类期刊数量:\n")
print(table(cluster_df$Cluster))

# 补充3.1.3
# ---- 作者共被引：快速出图 + 聚类表 ----
if (!("CR_AU" %in% colnames(M_all))) {
  M_all <- metaTagExtraction(M_all, Field = "CR_AU", sep = ";")
}

# 构建网络
auth_cocitation <- biblioNetwork(M_all, 
                                 analysis = "co-citation", 
                                 network = "authors", 
                                 sep = ";")

# 直接用 networkPlot 快速出图（只显示前30个高频作者）
png("Figure_Auth_Co-citation_Overall.png", width = 1200, height = 900, res = 150)
networkPlot(auth_cocitation, n = 30, 
            Title = "Author Co-citation Network (Top 30)", 
            type = "fruchterman", 
            labelsize = 0.6,
            size = TRUE)
dev.off()
cat("✅ 已生成: Figure_Auth_Co-citation_Overall.png\n")

# ---- 聚类分析（只导出表格，不绘图） ----
library(igraph)
g_auth <- graph_from_adjacency_matrix(auth_cocitation, weighted = TRUE, mode = "undirected")
g_auth <- delete.vertices(g_auth, degree(g_auth) == 0)

cl_auth <- cluster_louvain(g_auth)
cat("聚类数量:", length(unique(membership(cl_auth))), "\n")

# 聚类成员表
auth_cluster_df <- data.frame(
  Cluster = membership(cl_auth),
  Author = V(g_auth)$name
) %>%
  group_by(Cluster) %>%
  mutate(Rank_in_Cluster = row_number()) %>%
  ungroup() %>%
  arrange(Cluster, Rank_in_Cluster)

write.csv(auth_cluster_df, "Table_Auth_Co-citation_Clusters.csv", row.names = FALSE)

# 查看各聚类前5位作者
auth_cluster_df %>%
  group_by(Cluster) %>%
  slice_head(n = 5) %>%
  print(n = 50)

cat("✅ 已生成: Table_Auth_Co-citation_Clusters.csv\n")
auth_cluster_df <- read.csv("Table_Auth_Co-citation_Clusters.csv")
auth_cluster_df %>%
  group_by(Cluster) %>%
  slice_head(n = 5) %>%
  print(n = 50)

# ---- 3.2.1 作者合作网络（总体） ----
cat("【补跑】作者合作网络聚类分析...\n")

# 构建总体合作网络
auth_collab <- biblioNetwork(M_all, 
                             analysis = "collaboration", 
                             network = "authors", 
                             sep = ";")

# 转为igraph
library(igraph)
g_auth_collab <- graph_from_adjacency_matrix(auth_collab, weighted = TRUE, mode = "undirected")

# 剔除孤立节点
g_auth_collab <- delete.vertices(g_auth_collab, degree(g_auth_collab) == 0)

# 网络统计
cat("网络节点数:", vcount(g_auth_collab), "\n")
cat("网络边数:", ecount(g_auth_collab), "\n")
cat("网络密度:", edge_density(g_auth_collab), "\n")

# Louvain聚类
cl_collab <- cluster_louvain(g_auth_collab)
cat("聚类数量:", length(unique(membership(cl_collab))), "\n")

# 生成聚类成员表
collab_cluster_df <- data.frame(
  Cluster = membership(cl_collab),
  Author = V(g_auth_collab)$name,
  Degree = degree(g_auth_collab)
) %>%
  group_by(Cluster) %>%
  mutate(
    Cluster_Size = n(),
    Rank_in_Cluster = row_number(-Degree)  # 按度数排序
  ) %>%
  ungroup() %>%
  arrange(Cluster, Rank_in_Cluster)

# 查看各聚类规模
cluster_sizes <- collab_cluster_df %>%
  group_by(Cluster) %>%
  summarise(Size = n(), .groups = "drop") %>%
  arrange(desc(Size))

cat("\n各聚类规模:\n")
print(cluster_sizes)

# 查看前5大聚类的前10位成员
top_clusters <- cluster_sizes %>% slice_head(n = 5) %>% pull(Cluster)

cat("\n前5大聚类的核心成员:\n")
collab_cluster_df %>%
  filter(Cluster %in% top_clusters) %>%
  group_by(Cluster) %>%
  slice_head(n = 10) %>%
  print(n = 50)

# 保存完整聚类表
write.csv(collab_cluster_df, "Table_Collab_Network_Clusters.csv", row.names = FALSE)

# ---- 绘制合作网络图（前30个高频作者） ----
png("Figure_Auth_Collaboration_Overall.png", width = 1200, height = 900, res = 150)
networkPlot(auth_collab, n = 20, 
            Title = "Author Collaboration Network (Top 20)", 
            type = "fruchterman", 
            labelsize = 0.7,
            size = TRUE)
dev.off()
cat("✅ 已生成: Figure_Auth_Collaboration_Overall.png\n")

# ---- 3.2.2 国家分布与对比 ----
cat("【补跑】国家分布统计...\n")

# 提取总体国家字段
if (!("AU_CO" %in% colnames(M_all))) {
  M_all <- metaTagExtraction(M_all, "AU_CO", sep = ";")
}

# ----  总体国家排名 ----
country_all_raw <- strsplit(M_all$AU_CO, ";")
country_all <- unlist(country_all_raw)
country_all <- trimws(country_all)
country_all <- country_all[nchar(country_all) > 0]

country_df_all <- as.data.frame(sort(table(country_all), decreasing = TRUE))
names(country_df_all) <- c("Country", "Freq")
country_df_all$Percentage <- round(country_df_all$Freq / sum(country_df_all$Freq) * 100, 2)

write.csv(head(country_df_all, 20), "Table_Countries_Overall_Top20.csv", row.names = FALSE)
cat("✅ 已保存: Table_Countries_Overall_Top20.csv\n")

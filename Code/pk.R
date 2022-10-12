library(readxl)
library(tidyverse)
library(NonCompart)
library(plotly)
library("PerformanceAnalytics")
pkdata <- read_excel("Data/PK_NCA.xlsx")
pkdata %>% select(Cohort) %>% distinct()

pkdata1 <- pkdata %>% 
  gather(substance, conc, 4:6) %>%   
  mutate(conc = as.numeric(conc)) %>%
  filter(!is.na(conc))



ind_plot24 <- pkdata1 %>% filter(time <= 24) %>% 
  ggplot(aes(x=time, y=conc, col=substance)) +
  geom_line(alpha=0.3) +
  geom_point(alpha=0.3) +
  facet_wrap(~Subject, ncol=3, scales="free") +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Time (h)", y = "Plasma Concentration",
       col = "Substance",
       title = "Individual Concentration-Time Curves") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, colour = NA)) +
  ggsci::scale_color_aaas()

ggplotly(ind_plot24)

ind_plot <- pkdata1 %>% 
  ggplot(aes(x=time, y=conc, col=substance)) +
  geom_line(alpha=0.3) +
  geom_point(alpha=0.3) +
  facet_wrap(~Subject, ncol=3, scales="free") +
  scale_x_continuous(breaks = seq(0,72,24)) +
  labs(x = "Time (h)", y = "Plasma Concentration (?)",
       col = "Substance",
       title = "Individual Concentration-Time Curves") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, colour = NA)) +
  ggsci::scale_color_aaas()

# ggsave("ind_plot24.png", width = 12, height = 18)
# ggsave("ind_plot.png", width = 12, height = 18)


mean_plot24 <- pkdata1 %>% 
  filter(time <= 24) %>% 
  group_by(substance, time) %>% 
  mutate(mean=mean(conc), sd=sd(conc)) %>% 
  ggplot(aes(x=time, y=mean, col=substance)) +
  geom_point(size=2, alpha=0.5) + 
  geom_line(size=0.5, alpha=0.5) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd)) +
  facet_wrap(~substance) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  labs(x = "Time (h)", y = "Plasma Concentration (?)", col = " ") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, colour = NA)) +
  ggsci::scale_color_aaas()

mean_plot <- pkdata1 %>% 
  group_by(substance, time) %>% 
  mutate(mean=mean(conc), sd=sd(conc)) %>% 
  ggplot(aes(x=time, y=mean, col=substance)) +
  geom_point(size=2, alpha=1) + 
  geom_line(size=0.5, alpha=1) +
  geom_errorbar(aes(ymin=mean, ymax=mean+sd), alpha = 0.3, width = 0.2) +
  facet_wrap(~substance, scales="free_x") +
  scale_x_continuous(breaks = seq(0,72,24)) +
  labs(x = "Time (h)", y = "Plasma Concentration (?)", col = " ") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, colour = NA)) +
  ggsci::scale_color_aaas()
mean_plot
mean_plot24_both <- cowplot::plot_grid(mean_plot, mean_plot + scale_y_log10(), labels = LETTERS[1:2], ncol=1)
mean_plot24_both


# ggsave("mean_plot24.png", width = 15, height = 8) 
# ggsave("mean_plot.png", width = 15, height = 8) 
# ggsave("mean_plot24_both.png", width = 15, height = 10) 


ncadata1 <- pkdata1 %>% filter(time >= 0) %>% as.data.frame()

time0 <- pkdata1 %>% select(Subject, Cohort, substance) %>% distinct() %>% mutate(time = 0) %>% mutate(conc=0)

ncadata2 <- pkdata1 %>% mutate(time = time+2) %>%
  arrange(Subject, substance, time) %>% as.data.frame() 

nca1 <- tblNCA(ncadata1, 
              key=c("substance", "Subject"), 
              colTime="time", 
              colConc="conc",
              R2ADJ=0.1) %>%
        select(substance, Subject, CMAX, TMAX, LAMZHL, AUCLST, AUCIFO)

nca1 %>%
  group_by(substance) %>%
  summarize_at(vars(CMAX, TMAX, LAMZHL, AUCLST, AUCIFO), mean)
nca1


my_theme <-
  list(
    # round large p-values to two places
    "pkgwide-fn:pvalue_fun" = function(x) style_pvalue(x, digits = 2),
    "pkgwide-fn:prependpvalue_fun" = function(x) style_pvalue(x, digits = 2, prepend_p = TRUE),
    # report median (IQR) and n (percent) as default stats in `tbl_summary()`
    "tbl_summary-str:continuous_stat" = "{mean} ± {sd}",
    "tbl_summary-str:categorical_stat" = "{n} ({p}%)"
  )
set_gtsummary_theme(my_theme)
library(gtsummary)
library(kableExtra)
cordata_lm
density_cor <- function(i) {
  cordata_lm %>%
    mutate(sex = ifelse(sex == "남성", "male", "female")) %>%
    ggplot(aes(x = {{i}}, col = sex, fill = sex)) +
    geom_density(alpha = 0.3) +
    facet_wrap(vars(substance), scale = "free") +
    theme_bw() +
    theme(strip.background = element_rect(fill = NA, colour = NA)) +
    scale_color_manual(values = c("red", "blue"), aesthetics = c("colour", "fill"))
}
library(tidyverse)
density_cor(LAMZHL)
nca1 %>%
  select(-Subject, -TMAX) %>%
  tbl_summary(by=substance) %>%
  modify_header(label = "**Parameter**") %>%
  bold_labels()
library(gtsummary)
?geom_density
nca2 <- tblNCA(ncadata2, 
               key=c("substance", "Subject"), 
               colTime="time", 
               colConc="conc",
               adm="Infusion",
               dur=2,
               R2ADJ=0.1)

iddose <- read_excel("Data/SNB101_Phase I_raw data_20221011 (1).xlsx", sheet="dose") 
snrn1 <- read_excel("Data/SNB101_Phase I_raw data_20221011.xlsx", sheet="PKRaw_Cycle 1")  
snrn <- t(snrn1[1:2, 2:17]) %>% as.data.frame() %>%  distinct()
snrn
colnames(snrn) <- c("Subject", "SID")
colnames(iddose) <- c("SID", "dose")
iddose
iddose1 <- iddose %>% full_join(snrn) %>% mutate(dose = as.numeric(dose))
iddose1[1, 3] <- "C1_01-001"

dose_pk <- pkdata1 %>%
  left_join(iddose1, by = "Subject") 




nca <- nca1 %>% select(1:2, AUCLST, CMAX, LAMZHL) %>% mutate(method = 1) %>% 
  rbind(nca2 %>% select(1:2, AUCLST, CMAX, LAMZHL) %>% mutate(method = 2)) %>%  
  gather(param, raw, AUCLST:LAMZHL) %>% 
  left_join(iddose1) %>% 
  mutate(dn = raw/dose) %>% 
  gather(rowdn, value, c(raw, dn)) 
nca_dn1 <- nca1 %>%
  left_join(iddose1) %>%
  mutate(CMAX_dn = CMAX/dose, AUCLST_dn = AUCLST/dose, AUCIFO_dn = AUCIFO/dose) %>%
  select(-CMAX, -TMAX, -AUCLST, -AUCIFO) %>%
  select(Subject, dose, everything()) %>%
  rename(Dose = dose)

nca_dn1


pktable <- nca %>% 
  group_by(method, substance, param, rowdn) %>% 
  summarise_at(vars(value), lst(mean, sd, median, min, max)) %>% 
  ungroup() %>% 
  mutate_at(-(1:4), round, 2) %>% 
  arrange(method, rowdn, substance, param)


write_csv(pktable, "pktable.csv")  


dm <- read_excel("Data/SNB101_Phase I_raw data_20221011.xlsx", sheet="DM") %>% 
  select(1,4,5) 
hw <- read_excel("Data/SNB101_Phase I_raw data_20221011.xlsx", sheet = "HW") %>%
  select(1, 3, 6, 8)

colnames(dm) <- c("SID", "age", "sex")
colnames(hw) <- c("SID", "VISIT", "WT", "BSA")
cordata <- nca_dn1 %>% left_join(dm) %>% left_join(hw %>% filter(VISIT == "C1D1"))

colnames(cordata)
str(cordata)
cordata_r <- cordata %>%
  mutate(sex = ifelse(sex == "남성", 1, 0)) %>%
  mutate_at(vars("age", "WT", "BSA"), as.numeric) %>%
  select(substance, LAMZHL, CMAX_dn, AUCLST_dn, AUCIFO_dn, age, sex, WT, BSA) 

str(cordata_r)
chart.Correlation(cordata_r %>% filter(substance == "irinotecan") %>% select(-substance), histogram = TRUE, pch = 19)



cordata_lm <- cordata %>%
  mutate_at(vars("age", "WT", "BSA"), as.numeric) %>%
  select(substance, LAMZHL, CMAX_dn, AUCLST_dn, AUCIFO_dn, age, sex, WT, BSA)
m1 <- lm(LAMZHL ~ sex + age + WT+ BSA, data = cordata_lm)
tbl_regression(m1)
plot(cordata %>% filter(rowdn == "dn" & method == 1) %>% select(AUCLST, age, sex))
plot(cordata %>% filter(rowdn == "dn" & method == 1) %>% select(CMAX, age, sex))
plot(cordata %>% filter(rowdn == "dn" & method == 1 ) %>% select(LAMZHL, age, sex))

plot(cordata %>% filter(rowdn == "dn" & method == 2) %>% select(AUCLST, age, sex))
plot(cordata %>% filter(rowdn == "dn" & method == 2) %>% select(CMAX, age, sex))
plot(cordata %>% filter(rowdn == "dn" & method == 2 ) %>% select(LAMZHL, age, sex))

dev.off()


library(gtsummary)
library(tidyverse)

tbl <-
  c("cyl", "cyl + disp") %>% # vector of covariates
  map(
    ~ paste("mpg", .x, sep = " ~ ") %>% # build character formula
      as.formula() %>% # convert to proper formula
      lm(data = mtcars) %>% # build linear regression model
      tbl_regression() # display table with gtsummary
  ) %>%
  # merge tables into single table
  tbl_merge(
    tab_spanner = c("**Univariate**", "**Multivariable**")
  )
tbl
list(1, 2)

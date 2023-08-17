#### LOAD PACKAGES ####
pacman::p_load(tidyverse,readxl,gtsummary,corrplot,flextable,nlme,mgcv, 
               DHARMa,performance,emmeans,itsadug,gratia)

#### LOAD DATASET ####
DF = readxl::read_excel("COVID_PAMI_mod.xlsx") %>%
  mutate(
    # Tiempo a factor
    Tiempo_cat = factor(Tiempo),
    # Edad a factor
    Edad_cat = gtools::quantcut(Edad, q = 2),
    # UA logaritmica
    UA_log = log(UA_mod),
    # Tratamiento (grupo:vacuna)
    Tratamiento = fct_cross(Grupo,Vacuna) %>% 
      fct_relevel("Naïve:Sinopharm","Naïve:Sputnik V", after = 1),
    # Variables categóricas a factor
    across(where(is.character), as.factor)) %>% 
  # Elimina registros con datos faltantes para comorbilidades
  filter(!is.na(Diabetes)) %>% droplevels()

### Descripción grupos----
DF %>% filter(Tiempo==0) %>% select(Grupo, Vacuna) %>% 
  tbl_cross() %>% bold_labels()

### Descripción edades----
compareGroups::descrTable(Sexo ~ Edad_cat, data = DF, method = 4, show.all = T)

### Spaghetti plots----
ggpubr::ggarrange(
  # variable respuesta: UA
  ggplot(DF, aes(x = factor(Tiempo), y = UA_mod, group = ID, color = Tratamiento)) +
    geom_line() +
    labs(y = "UA", x = "Tiempo (días)") +
    scale_color_grey(start = .4, end = 0) +
    # scale_color_brewer(palette = "Set2") +
    theme_light()
  ,
  # Variable respuesta: log(UA)
  ggplot(DF, aes(x = factor(Tiempo), y = UA_log, group = ID, color = Tratamiento)) +
    geom_line() +
    labs(y = "log(UA)", x = "Tiempo (días)") +
    scale_color_grey(start = .4, end = 0) +
    # scale_color_brewer(palette = "Set2") +
    theme_light(), 
  common.legend = T, legend = "bottom", labels = "AUTO",font.label = list(size=12))

### Interacción entre grupo y vacuna----
fit = lme(UA_log ~ Grupo*Vacuna, random = ~1|ID, data = DF)
anova(fit)

### FIG. 1: Gráfico de perfiles----
g1 = ggplot(data = DF, aes(x = Tiempo, y = UA_log, color = Tratamiento, group = Tratamiento)) +
  geom_smooth(method = "loess", se = F) + 
  # Colorblind friendly
  scale_color_brewer(palette = "Set2") +
  # # Greyscale
  # scale_color_grey(end = .9) +
  labs(y = "log(UA)") +
  scale_x_continuous(breaks = c(0,21,42,120,180), minor_breaks = F) +
  theme_light() + theme(legend.title = element_blank())

# Gráfico con leyenda
g1 + theme(legend.position = "bottom")

## Guarda gráfico
# ggsave("FIGS/Fig1.svg", width = 16, height = 10, units = "cm", dpi = 300)
# ggsave("FIGS/Fig1bw.svg", width = 16, height = 10, units = "cm", dpi = 300)

# # Gráfico sin leyenda
# g1 + theme(legend.position = "none")
# 
# ## Guarda gráfico
# ggsave("FIGS/Fig1sl.svg", width = 16, height = 10, units = "cm", dpi = 300)
# 
# # plot legend
# require(grid)
# grid.newpage()
# grid.draw(cowplot::get_legend(g1))

### Tabla 1: Asociación con VE----
DF %>% select(Sexo, Edad_cat, COVID_estudio, Brote_estudio, Diabetes, 
              HTA, EPOC, Obesidad_severa, Enf_renal_cronica, Insuf_cardiaca,
              Inmunodeficiencia, UA_log, ID) %>% 
  tbl_uvregression(y = UA_log, method = glmmTMB::glmmTMB,
                   formula = "{y} ~ {x} + (1|ID)",
                   pvalue_fun = function(x) style_pvalue(x, digits = 2),
                   show_single_row = -UA_log,
                   label = list(Sexo ~ "Sexo: Masculino",
                                Edad_cat ~ "Grupo etario",
                                COVID_estudio ~ "COVID durante estudio: Si",
                                Brote_estudio ~ "Brote durante el estudio: Si",
                                Diabetes ~ "Diabetes: Si",
                                HTA ~ "HTA: Si",
                                EPOC ~ "EPOC: Si",
                                Obesidad_severa ~ "Obesidad severa: Si",
                                Enf_renal_cronica ~ "Enf, renal crónica: Si",
                                Insuf_cardiaca ~ "Insuf, cardiaca: Si",
                                Inmunodeficiencia ~ "Inmunodeficiencia: Si")) %>% 
  bold_p() %>% 
  as_flex_table() %>% 
  set_caption("Tabla 1. Resultados de los modelos lineales mixtos (LMM) univariados para selección de variables explicativas. Beta: coeficiente de regresión; 95% IC: intervalo de confianza al 95%. Los P-valores estadísticamente significativos se muestran en negritas.") %>% 
  save_as_docx(path = "FIGS/tab1.docx")

### Asociación entre VE----
require(compareGroups)

# Sexo
descrTable(Sexo ~ Edad_cat + # asoc. significativa
             COVID_estudio +
             Brote_estudio +
             Enf_renal_cronica + # asoc. significativa
             Insuf_cardiaca + 
             Inmunodeficiencia, 
           data = DF %>% filter(Tiempo==0))

# Grupo etario
descrTable(Edad_cat ~ COVID_estudio +
             Brote_estudio +
             Enf_renal_cronica +
             Insuf_cardiaca + 
             Inmunodeficiencia, 
           data = DF %>% filter(Tiempo==0))

# COVID durante el estudio
descrTable(COVID_estudio ~ Brote_estudio + # asoc. significativa
             Enf_renal_cronica + 
             Insuf_cardiaca + # asoc. significativa
             Inmunodeficiencia, 
           data = DF %>% filter(Tiempo==0))

# Brote durante el estudio
descrTable(Brote_estudio ~ Enf_renal_cronica + # asoc. significativa
             Insuf_cardiaca + 
             Inmunodeficiencia, 
           data = DF %>% filter(Tiempo==0))

# Enfermedad renal crónica
descrTable(Enf_renal_cronica ~ Insuf_cardiaca + # asoc. significativa
             Inmunodeficiencia, 
           data = DF %>% filter(Tiempo==0))

# Insuficiencia cardíaca
descrTable(Insuf_cardiaca ~ Inmunodeficiencia, 
           data = DF %>% filter(Tiempo==0))

### FIG. S2: Autocorrelación temporal----
require(corrplot)
# svg(filename = "FIGS/FigS2.svg", width = 4, height = 4)
# svg(filename = "FIGS/FigS2bw.svg", width = 4, height = 4)
DF %>% select(ID, Tiempo, UA_log) %>% 
  pivot_wider(values_from = "UA_log", names_from = "Tiempo") %>% select(-ID) %>% 
  cor(., use = "na.or.complete") %>% 
  corrplot(type = "upper", addCoef.col = "black",outline = F,
           diag = F, 
           mar = c(1,1,1,1), tl.col = "black",
           # Colorblind-friendly
           # col = COL2("PuOr"))
           # Greyscale
           col = COL1("Greys"))
# dev.off()

#### MODELOS MARGINALES (GLS) ####
### Autorregresiva continua de orden 1 c/varianza homogénea
gls1 = gls(UA_log ~ Tratamiento*Tiempo_cat +
             Edad_cat + Brote_estudio + Insuf_cardiaca + Inmunodeficiencia,
           correlation = corCAR1(form = ~1|ID), 
           data = DF)

### Autorregresiva continua de orden 1 c/varianza heterogénea
gls2 = gls(UA_log ~ Tratamiento*Tiempo_cat +
             Edad_cat + Brote_estudio + Insuf_cardiaca + Inmunodeficiencia,
           correlation = corCAR1(form = ~1|ID), 
           weights = varIdent(form= ~ 1|Tiempo_cat),
           data = DF)

### Compara modelos marginales----
compare_performance(gls1, gls2, metrics = "all", rank = F)

### Selección variables explicativas----
anova(gls2)
gls2a = update(gls2, ~.-Insuf_cardiaca)
anova(gls2a)
gls2b = update(gls2a, ~.-Edad_cat)
anova(gls2b)
gls2c = update(gls2b, ~.-Inmunodeficiencia)
anova(gls2c)

### Tabla 2: Selección de variables
MuMIn::model.sel(gls2, gls2a, gls2b, gls2c, rank = "AIC") %>% 
  as_tibble() %>% 
  mutate(Modelo = c("- Inmunodeficiencia","- Grupo etario", "- Insuf. cardíaca", "Saturado")) %>% 
  select(Modelo, df,logLik, AIC, delta) %>%
  flextable(cwidth = c(1.5,.5,1,1,1)) %>% 
  set_caption("Tabla 2. Selección de variables explicativas en el modelo GLS con estructura de correlación temporal autorregresiva continua de primer orden (corCAR1) y heterogeneidad de varianzas (varIdent).") %>% 
  save_as_docx(path = "./FIGS/tab2.docx")

### Limpia environment
rm(list = setdiff(ls(),c("DF","gls2c")))

### Análisis residuales----
check_model(gls2c)

### Test de comparaciones múltiples----
em = emmeans(gls2c, ~ Tratamiento * Tiempo_cat + Brote_estudio)

# Resumen
pairs(em, infer = T, simple = "Tratamiento") %>% 
  as.data.frame() %>% 
  mutate_if(is.numeric, round, 3) %>% 
  flextable() %>% save_as_html(path = "emmeans_tab.html")

## FIG. S3: Presentación resultados test comparaciones múltiples----
g1 = ggeffects::ggemmeans(model = gls2c,
                          terms = c("Tiempo_cat","Tratamiento","Brote_estudio"))

ggplot(g1, aes(x = x %>% as.character %>% as.numeric, y = predicted, color = group)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_line(aes(lty = group)) +
  # Colorblind friendly
  # scale_color_brewer(palette = "Set2") +
  # Greyscale
  scale_color_grey(start = .9, end = .2) +
  facet_wrap(.~ facet, ncol = 1,
             labeller = labeller(facet = c(No = "Brote estudio: No", Si = "Brote estudio: Si"))) +
  scale_x_continuous(breaks = c(0,21,42,120,180)) +
  labs(x = "Tiempo (días)", y = "log(UA)") +
  theme_minimal() + theme(legend.position = "bottom", legend.title = element_blank(),
                          axis.title = element_text(face = "bold", size = 10),
                          strip.text = element_text(face = "bold", size = 10))

# ggsave("FIGS/FigS3.svg", width = 16, height = 20, units = "cm", dpi = 300)
# ggsave("FIGS/FigS3bw.svg", width = 16, height = 20, units = "cm", dpi = 300)

#### MODELA RESPUESTA NO LINEAL (GAMM) ####
### GAMM c/corCAR1----
gamm = gamm(UA_log ~ Brote_estudio + Tratamiento +
              s(Tiempo, by = Tratamiento, k = 5),
            random = list(ID = ~ 1),
            correlation = corCAR1(form = ~1|ID), 
            weights = varIdent(form= ~ 1|Tiempo_cat),
            data = DF, method = "REML")

summary(gamm$gam)
summary(gamm$lme)

### Análisis residuales----
acf(gamm$lme %>% residuals())
pacf(gamm$lme %>% residuals())

appraise(gamm$gam)

### Grados de libertad efectivos----
edf(gamm$gam)

### FIG. 2: Grafica curvas suavizadas----
require(gridExtra)
svg(filename = "FIGS/Fig2sl.svg", width = 6.3, height = 7.9)
par(mfrow = c(2,1), cex = .75)
# Brote: Si
plot_smooth(gamm$gam, view = "Tiempo", plot_all = "Tratamiento", 
            cond = list(Brote_estudio="Si"), v0 = c(0,21,42,120,180),
            rug = F, rm.ranef = F, legend_plot_all = F,
            ylim = c(2,11), xlim = c(0,200), 
            ylab = "log(UA)", hide.label = F, main = "(A)", adj = 0,
            # Colorblind friendly
            col = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F"))
            

# Brote: No
plot_smooth(gamm$gam, view = "Tiempo", plot_all = "Tratamiento", 
            cond = list(Brote_estudio="No"), v0 = c(0,21,42,120,180),
            rug = F, rm.ranef = F,legend_plot_all = F,
            ylim = c(2,11), xlim = c(0,200),
            ylab = "log(UA)", hide.label = F, main = "(B)", adj = 0,
            # Colorblind friendly
            col = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F"))
dev.off()

### FIG. 3: Gráfico de los términos paramétricos----
# svg(filename = "FIGS/Fig3.svg", width = 6.3, height = 7.9)
par(mfrow = c(2,1), cex = .75)
plot_parametric(gamm$gam, pred = list(Tratamiento = levels(DF$Tratamiento)),
                cond = list(Brote_estudio="Si"), xlim = c(4, 12),
                xlab = "log(UA)", main = "")
title(main = "(A)", adj = 0)

plot_parametric(gamm$gam, pred = list(Tratamiento = levels(DF$Tratamiento)),
                cond = list(Brote_estudio="No"),  xlim = c(4, 12),
                xlab = "log(UA)", main = "")
title(main = "(B)", adj = 0)


# dev.off()


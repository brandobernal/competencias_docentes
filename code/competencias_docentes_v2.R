# paqueterias -------------------------------------------------------------

rm(list=ls())


options(scipen = 999)

pacman::p_load(tidyverse, psych, lavaan, semTools, semPlot, here, MVN, nortest, mirt)

# bases de datos ----------------------------------------------------------

bd <- read.csv(here("datos", "bd_competencias_docentes.csv"))
bd <- bd %>% rename(sexo = ï..Sexo)
bd %>% names()

# calculo de poder muestral

power_size <- function(n, z_score, proportion, error){
  pwr_dw_up <- z_score^2*(proportion*(1-proportion))
  pwr_dw_dw <- error^2*(n)
  pwr_dw = 1+(pwr_dw_up/pwr_dw_dw)
  pwr_up = pwr_dw_up/error^2
  a<-round(pwr_up/pwr_dw, 0)
  message("\nLa muestra mínima para ", {n}, " casos es de:"); print(a)
  message("con un error de estimación de +- 5 casos")}

power_size(27227, 1.96, .5, .05)

# tratamiento de los datos ------------------------------------------------

# reconfigurar las puntuaciones

bd[2:18] %>% lapply(table) # hay un 21 en $A2
bd[bd == 21] <- 2 #reconfigurar el valor

# items 12, 13 y 14 se califican al reves
head(bd[12:14])
cols = c("A1", "A2", "A3"); bd[ ,cols] = 6 - bd[ ,cols]; head(bd[12:14])

# obtener los indices de ambas escalas

bd$ig_eecd <- bd[2:10] %>% rowSums()
bd$ig_iscs <- bd[11:18] %>% rowSums() 

bd[19:20] %>% lapply(table)

# asignar nivel de medicion de las variables

bd %>% sapply(class)
bd$sexo <- as_factor(bd$sexo); bd %>% sapply(class)

# buscar outliers

boxplot(bd$ig_eecd, 
        main = "Indice general EECD de los estudiantes",
        xlab = "Estudiantes",
        ylab = "Valor del indice general"); boxplot.stats(bd$ig_eecd)

bd <- bd %>% filter(ig_eecd > 21) # eliminar outliers 

boxplot(bd$ig_iscs, 
        main = "Indice general ISCS de los estudiantes",
        xlab = "Estudiantes",
        ylab = "Valor del indice general"); boxplot.stats(bd$ig_iscs)

# buscar valores perdidos

bd %>% sapply(function(x) sum(is.na(x)));
bd %>% apply(2, (function(x){sum(is.na(x))/length(x)*100}))

# analisis descriptivos ---------------------------------------------------

hist(bd$ig_eecd, breaks = 15, 
     main = "Histograma para EECD", xlab = "Valores de IG", ylab = "Frecuencia")

hist(bd$ig_iscs, breaks = 15, 
     main = "Histograma para ISCS", xlab = "Valores de IG", ylab = "Frecuencia")

# estadisticos de la muestra ----------------------------------------------

bd %>% 
  group_by(sexo) %>% 
  summarise(n = n(),
            IG_eecd = mean(ig_eecd),
            DE_eecd = sd(ig_eecd),
            IG_iscs = mean(ig_iscs),
            DE_iscs = sd(ig_iscs))

# graficos de densidad para ig de eecd y iscs

ggplot(na.omit(bd), aes(x = ig_eecd)) +
  geom_density(aes(group = sexo, fill = sexo), alpha = 0.4) + 
  scale_fill_discrete(name = "Grupos", 
                      labels = c("Mujer", "Hombre"))

ggplot(na.omit(bd), aes(x = ig_iscs)) +
  geom_density(aes(group = sexo, fill = sexo), alpha = 0.4) + 
  scale_fill_discrete(name = "Grupos", 
                      labels = c("Mujer", "Hombre"))

# boxplot en funcion de sexo para ig de eecd y iscs

bd %>% ggplot() + 
  geom_boxplot(aes(x = sexo, y = ig_eecd, group=sexo, fill = sexo)) +
  labs(x = "Sexo", y = "Valores de IG", 
       title = "Puntuaciones de IG EECD en función a sexo", 
       subtitle = "Separados por sexo") +
  scale_fill_hue(labels = c("Mujer", "Hombre"))

bd %>% ggplot() + 
  geom_boxplot(aes(x = sexo, y = ig_iscs, group=sexo, fill = sexo)) +
  labs(x = "Sexo", y = "Valores de IG", 
       title = "Puntuaciones de IG ISCS en función a sexo", 
       subtitle = "Separados por sexo") +
  scale_fill_hue(labels = c("Mujer", "Hombre"))

# estadisticos preliminares -----------------------------------------------

# eecd

bd[2:10] %>% describe()
bd[2:10] %>% alpha() # cronbach
bd[2:10] %>% omega(nfactors = TRUE) # omega mcdonald

# iscs

bd[11:18] %>% describe()
bd[11:18] %>% alpha() # cronbach
bd[11:18] %>% omega(nfactors = TRUE) # omega mcdonald

# correlaciones -----------------------------------------------------------

# eecd

bd[2:10] %>% cor() %>% # correlaciones bivariadas
    knitr::kable(digits = 2, caption = "Tabla N. Correlaciones entre items de la escala")

bd[2:10] %>% cor() %>% KMO() # kmo

# iscs

bd[11:18] %>% cor() %>% # correlaciones bivariadas
  knitr::kable(digits = 2, caption = "Tabla N. Correlaciones entre items de la escala")

bd[11:18] %>% cor() %>% KMO() # kmo


# normalidad --------------------------------------------------------------

# eecd

lillie.test(bd$ig_eecd) # p-value < 0.05, no se asume normalidad

bd[2:10] %>% mvn(mvnTest = "hz") # no normalidad multivariante; estimador dwls robusto

# iscs

lillie.test(bd$ig_iscs) # p-value < 0.05, no se asume normalidad

bd[11:18] %>% mvn(mvnTest = "hz") # no normalidad multivariante; estimador dwls robusto

# dimensionalidad ---------------------------------------------------------

CFA <- function(datos, modelo, estimador){
  fit <- cfa(model = modelo, data = datos, estimator = "WLSMV")
  a <- fit %>% parameterEstimates(standardized = TRUE, se = TRUE, zstat = FALSE, 
                                  pvalue = FALSE, ci = FALSE)
  a1 <- a[, -c(4:5, 8)] %>% knitr::kable(.,
                                         col.names = c("factor" ,"", "item", "latent std", "all std"))
  b <- fit %>% fitmeasures(c("chisq", "df", "pvalue", "cfi", "tli", 
                             "gfi", "nfi", "rmsea", "rmsea.ci.lower", 
                             "rmsea.ci.upper", "srmr","aic", "bic"))
  b1 <- b %>% knitr::kable(digits = options(digits = 2),
                           col.names = c("valor"),
                           caption = "indices de ajuste para el modelo")
  c <- fit %>% semPaths("par", "std", weighted = FALSE, nCharNodes = 7, 
                        shapeMan = "rectangle", sizeMan = 5, sizeMan2 = 2, 
                        rotation = 2)
  options(warn=-1)
  message("\nLos datos del modelo son: \n"); print(a)
  message("\nLos indices de ajuste del modelo son: \n"); print(b)}

# cfa ---------------------------------------------------------------------

# eecd

eecd_uni_fac <- {"# measurement model
            f1 =~ P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9"}

  # mod_eecd <- eecd_uni_fac %>% cfa(bd, "WLSMV")

mod_eecd <- readRDS(here("results", "mod_eecd.RDS"))
mod_eecd %>% lavInspect(what='rsquare')
mod_eecd %>% fitmeasures()

CFA(bd, eecd_uni_fac)

  # los resultados no coinciden con el modelo de la bd antigua

mod_antiguo_eecd <- readRDS(here("results", "mod_antiguos", "mod_antiguo_eecd.rds"))

fitMeasures(mod_antiguo_eecd)

mod_antiguo_eecd %>% parameterEstimates(standardized = TRUE, se = TRUE, zstat = FALSE, 
                   pvalue = FALSE, ci = FALSE)

#mod_antiguo_eecd %>% fitmeasures(c("chisq", "df", "pvalue", "cfi", "tli", 
 #             "gfi", "nfi", "rmsea", "rmsea.ci.lower", 
 #             "rmsea.ci.upper", "srmr","aic", "bic"))

comparacion_mods_eecd <- readRDS(here("results", "comparacion_mods_eecd.rds"))

# iscs

iscs_dos_fac <- {"# measurement model
            f1 =~ S1 + S5 + S6 + S7 + S8
            f2 =~ A1 + A2 + A3"}

  # mod_iscs <- cfa(iscs_dos_fac, bd, "WLSMV")

readRDS(here("results", "mod_iscs.RDS"))

iscs_model <- bd %>% CFA(iscs_dos_fac)

# invarianza --------------------------------------------------------------

invarianza <- function(datos = bd, modelo, variable){
  config <- cfa(model = modelo, data = na.omit(datos), group = variable)
  debil <- cfa(model = modelo, data = na.omit(datos), group = variable, 
               group.equal="loadings")
  fuerte <- cfa(model = modelo, data = na.omit(datos), group = variable, 
                group.equal=c("loadings", "intercepts"))
  estricto <- cfa(model = modelo, data = na.omit(datos), group = variable, 
                  group.equal=c("loadings", "intercepts", "residuals"))
  inv <- compareFit(config, debil, fuerte, estricto)
  options(warn=-1)
  message("\nEstadísticos de invarianza para la variable ", variable,"\n"); print(inv)}

# invarianza sexo para eecd

inv_sexo_eecd <- invarianza(bd, eecd_uni_fac, "sexo")

inv_sexo_eecd  <- readRDS(here("results", "inv_sexo_eecd.rds"))

  # no coincide con la invarianza de v1

inv_sexo_antiguo <- readRDS(here("results", "inv_sexo_antigua.rds"))

list(antigua = inv_sexo_antiguo,
     actual = inv_sexo)

# invarianza sexo para iscs

inv_sexo_iscs <- invarianza(bd, iscs_dos_fac, "sexo")

inv_sexo_iscs <- readRDS(here("results", "inv_sexo_iscs.rds"))

# modelos de la tri -------------------------------------------------------

# modelos eecd

items <- ifelse(bd[2:10]>= 5, 1, 0); head(items, 1) #reconfigurar respuestas

M2_mod_1PL <- M2(mirt(items, 1, itemtype = "Rasch"))
M2_mod_2PL <- M2(mirt(items, 1, itemtype = "2PL"))
M2_mod_3PL <- M2(mirt(items, 1, itemtype = "3PL"))

comparacion_tri <-("1PL" = M2_mod_1PL) %>% 
                      bind_rows("2PL" = M2_mod_2PL) %>% 
                      bind_rows("3PL" = M2_mod_3PL)

comparacion_tri <- readRDS(here("results", "comparacion_tri.rds"))

  # no coincide con los modelos tri antiguos de la eecd
  
comparacion_tri_antiguo <- readRDS(here("results", "mods_antiguos", "comparacion_tri_antiguo.rds"))

list(antiguo = comparacion_tri_antiguo, actual = comparacion_tri)

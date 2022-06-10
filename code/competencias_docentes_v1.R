# paqueterias -------------------------------------------------------------

rm(list=ls())

pacman::p_load(tidyverse, psych, lavaan, semTools, semPlot, here, MVN, mirt)

# bases de datos ----------------------------------------------------------

bd <- read.csv(here("datos", "datos_antiguos.csv"))
bd$ig <- bd[2:10] %>% rowSums()
bd <- bd %>% rename(sexo = ï..Sexo)

# asignar nivel de medicion de las variables

bd %>% sapply(class)
bd$sexo <- as_factor(bd$sexo); bd %>% sapply(class)

# calculo de poder muestral

power_size <- function(n, z_score, proportion, error){
  pwr_dw_up <- z_score^2*(proportion*(1-proportion))
  pwr_dw_dw <- error^2*(n)
  pwr_dw = 1+(pwr_dw_up/pwr_dw_dw)
  pwr_up = pwr_dw_up/error^2
  a<-round(pwr_up/pwr_dw, 0)
  message("\nLa muestra mínima para ", {n}, " casos es de:"); print(a)
  message("con un error de medicion de +- 5 casos")}

power_size(27227, 1.96, .5, .05)

# tratamiento de los datos ------------------------------------------------

# no hay outliers

boxplot(bd$ig, 
        main = "Indice general de los estudiantes",
        xlab = "Estudiantes",
        ylab = "Valor del indice general"); boxplot.stats(bd$ig)

# no hay valores perdidos

bd %>% sapply(function(x) sum(is.na(x)));
bd %>% apply(2, (function(x){sum(is.na(x))/length(x)*100}))

# analisis descriptivos ---------------------------------------------------

bd[2:10] %>% lapply(table) # puntuaciones de los items

hist(bd$ig, breaks = 10, 
     main = "Histograma para EECD", xlab = "Valores de IG", ylab = "Frecuencia")

# estadisticos de la muestra ----------------------------------------------

bd %>% 
  group_by(sexo) %>% 
  summarise(n = n(),
            IG = mean(ig),
            DE = sd(ig),
          )

ggplot(na.omit(bd), aes(x = ig)) +
  geom_density(aes(group = sexo, fill = sexo), alpha = 0.4) + 
  scale_fill_discrete(name = "Grupos", 
                      labels = c("Mujer", "Hombre"))
bd %>% ggplot() + 
  geom_boxplot(aes(x = sexo, y = ig, group=sexo, fill = sexo)) +
  labs(x = "Sexo", y = "Valores de IG", 
       title = "Puntuaciones de IG en funcion a sexo", 
       subtitle = "Separados por sexo") +
  scale_fill_hue(labels = c("Mujer", "Hombre"))

# estadisticos preliminares -----------------------------------------------

bd[2:10] %>% describe()

bd[2:10] %>% alpha() # cronbach

bd[2:10] %>% omega(nfactors = TRUE) # omega mcdonald

# correlaciones

bd[2:10] %>% 
  cor() %>% 
    knitr::kable(digits = 2, caption = "Tabla N. Correlaciones entre items de la escala")

bd[2:10] %>% cor() %>% KMO()

bd[2:10] %>% mvn(mvnTest = "hz") # no normalidad multivariante; estimador dwls robusto

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
  d <- fit %>% lavInspect(what='rsquare')
  options(warn=-1)
  message("\nLos datos del modelo son: \n"); print(a)
  message("\nLos indices de ajuste del modelo son: \n"); print(b)
  message("\nLa varianza explicada para cada item es: \n"); print(d)}

mod_uni_fac <- {"# measurement model
            f1 =~ P21EF + P22EF + P23EF + P24EF + P25EF + P26EF + P27EF + P28EF + P29EF"}

uni_fac <- bd %>% CFA(mod_uni_fac)

uni_fac <- readRDS(here("results", "uni_fac.rds"))

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
  message("\nEstadísticos de invarianza para la variable ", variable,"\n"); return(inv)}

inv_sexo <- invarianza(bd, mod_uni_fac, "sexo")

inv_sexo <- readRDS(here("results", "inv_sexo.rds"))

# modelos de la TRI -------------------------------------------------------

table(bd$P21EF)

items <- ifelse(bd[2:10]>= 5, 1, 0); head(items)

mod_1PL<- mirt(items, 1, itemtype = "Rasch")

M2_mod_1PL <- M2(mirt(items, 1, itemtype = "Rasch"))
M2_mod_2PL <- M2(mirt(items, 1, itemtype = "2PL"))
M2_mod_3PL <- M2(mirt(items, 1, itemtype = "3PL"))

a <- M2_mod_1PL %>% 
  bind_rows(M2_mod_2PL) %>% 
  bind_rows(M2_mod_3PL)
  
comparacion_tri.rds <- readRDS(here("results", "comparacion_tri.rds")); comparacion_tri.rds

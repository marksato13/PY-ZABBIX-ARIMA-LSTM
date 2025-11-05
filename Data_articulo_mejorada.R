attach(tesita_zabbix)
names(tesita_zabbix)
install.packages("tseries")
install.packages("astsa")
install.packages("forecast")
install.packages("foreign")
install.packages("quantmod")

library(astsa)
library(tseries)
library(lubridate)
library(tidyverse)
library(forecast)

cpu <- tesita_zabbix$cpu_usage
cpu_ts <- ts(cpu, start = c(1), frequency = 720)
cpu_ts
plot(cpu_ts)
plot(cpu_ts, main="Serie temporal del uso de CPU", ylab="CPU Usage (%)", xlab="Día")
#Observamos que de lo que está en lo mínimo sube y luego baja a full, así que tenemos que hacer estacionariedad
#ya que es requisito para hacer el ModeloArima, lo siguiente que haremos es para estacionar mediante log
#No aplicaremos transformación logarítmica porque necesitamos los valores reales del CPU para la predicción

#Limpiamos outliers sin distorcionar la escala, luego se elimina la tendencia
cpu_clean <- tsclean(cpu_ts)        # limpia outliers sin distorsionar escala
cpu_diff  <- diff(cpu_clean, 1)     # elimina tendencia
#Prueba de ADF para ver si es estacionaria
adf.test(cpu_clean)
adf.test(cpu_diff)

#Ahora haremos la función de auto-correlación y la de auto-correlación parcial, sirve para 
#saber cuantas medias móviles y cuantos auto-regresivos vamos a utilizar en nuestro modelo ARIMA
par(mfrow=c(2,1), mar=c(4,4,4,1)+.1)
acf(cpu_diff, lag.max = 100, main="ACF del uso de CPU (diferenciada)")
pacf(cpu_diff, lag.max = 100, main="PACF del uso de CPU (diferenciada)")
#Por si quieres restaurar para que no aparezca la gráfica 2 en 1 en los siguientes codes:
par(mfrow=c(1,1))


# Ajuste con p=1, d=1, q=0 sobre la serie limpia (no la diferenciada)
modelo1 <- arima(cpu_clean, order = c(1,1,0))
summary(modelo1)
tsdiag(modelo1)

modelo2 <- arima(cpu_clean, order = c(1,1,1))
summary(modelo2)
tsdiag(modelo2)

modelo3 <- arima(cpu_clean, order = c(2,1,0))
summary(modelo3)
tsdiag(modelo3)

modelo4 <- arima(cpu_clean, order = c(2,1,1))
summary(modelo4)
tsdiag(modelo4)

modelo5 <- arima(cpu_clean, order = c(3,1,3))
summary(modelo5)
tsdiag(modelo5)

modelo6 <- arima(cpu_clean, order = c(3,1,2))
summary(modelo6)
tsdiag(modelo6)
#
checkresiduals(modelo6)


#PRUEBAAAS
y <- cpu_clean
# Helper: ajusta, extrae métricas y diagnostics
fit_one <- function(order, lb_lag = 24) {
  m <- Arima(y, order = order, include.constant = FALSE)
  
  # AICc manual
  k <- length(coef(m))               # nº de parámetros estimados (p+q [+ estacionales si hubiera])
  n <- length(y)
  AIC_val  <- AIC(m)
  AICc_val <- AIC_val + (2 * k * (k + 1)) / (n - k - 1)
  
  # Varianza de los residuos (equivalente a sigma^2 del modelo)
  sig2 <- mean(residuals(m)^2, na.rm = TRUE)
  
  # Ljung–Box en lags bajos (importa para ruido blanco)
  lb <- Box.test(residuals(m), lag = lb_lag, type = "Ljung-Box", fitdf = k)
  
  data.frame(
    model  = paste0("ARIMA(", paste(order, collapse=","), ")"),
    AICc   = AICc_val,
    BIC    = BIC(m),
    sigma2 = sig2,
    LB_p   = lb$p.value,
    stringsAsFactors = FALSE
  )
}

# Candidatos que quieres probar
cands <- list(
  c(2,1,2), c(3,1,1), c(3,1,2), c(3,1,3),
  c(2,1,1), c(1,1,1), c(2,1,0), c(1,1,0)
)

tab <- do.call(rbind, lapply(cands, fit_one))
tab[order(tab$AICc), ]  # ordenado por AICc (menor es mejor)


#AutoArima para comparar
modelo7_auto <- auto.arima(
  y,
  seasonal       = FALSE,          # evita SARIMA con periodo 720
  stepwise       = FALSE,          # búsqueda global (mejor)
  approximation  = FALSE,          # sin aproximaciones
  ic             = "aicc",         # usamos AICc para comparar
  allowmean      = FALSE,          # consistente con tus ARIMA(d>=1) sin constante
  allowdrift     = TRUE            # permite drift si d=1 (lo decidirá el algoritmo)
)
summary(modelo7_auto)

# 3) Diagnóstico de residuos (mejor que tsdiag para objetos 'forecast')
checkresiduals(modelo7_auto)
# Probar con menos lags
Box.test(residuals(modelo7_auto), lag = 20, type = "Ljung-Box", fitdf = 5)

# Para ver la diferencia de la escala del MASE
  # 1. Predicciones dentro de la muestra
pred <- fitted(modelo7_auto)

  # 2. Calcular MAE del modelo
mae_model <- mean(abs(cpu_clean - pred))

  # 3. Calcular MAE del modelo naive (y_{t} - y_{t-1})
naive_pred <- cpu_clean[-length(cpu_clean)]
mae_naive <- mean(abs(diff(cpu_clean)))

  # 4. Calcular MASE manual
mase_manual <- mae_model / mae_naive
mase_manual


# Pronóstico de 1 día (720 muestras)
pronostico1 <- forecast(modelo6, h = 720)

# Ver valores y rangos
pronostico1

# Gráfico (forecast::autoplot si tienes ggplot2; si no, plot(fc_1d))
autoplot(pronostico1) +
  ggtitle("Pronóstico de uso de CPU (ARIMA(3,1,2)) – Próximo día") +
  xlab("Tiempo (min)") + ylab("CPU Usage (%)")


# Algo mas realista si queremos que varíe agregando tendencia o drift
modelo_drift <- Arima(cpu_clean, order = c(3,1,2), include.drift = TRUE)
fc_drift <- forecast(modelo_drift, h = 720)
autoplot(fc_drift) +
  ggtitle("Pronóstico con drift (ARIMA(3,1,2)+drift)") +
  ylab("CPU Usage (%)")
summary(modelo_drift)

# Índices/tiempos del horizonte futuro
t_future <- time(fc_drift$mean)

# --- Opción A: línea roja de umbral ---
autoplot(fc_drift) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Pronóstico con umbral 80% (línea)",
       y = "CPU Usage (%)", x = "Time")

# Comparando AIC
AIC(modelo6, modelo_drift)# ARIMA(3,1,2 y ARIMA(3,1,2)+drift


# Validando en un día hold-out (el 7mo día real)
n_per_day <- 720
y <- cpu_clean
train <- head(y, length(y)-n_per_day)
test  <- tail(y, n_per_day)

m_base   <- Arima(train, order=c(3,1,2))
m_drift  <- Arima(train, order=c(3,1,2), include.drift=TRUE)

fc_base  <- forecast(m_base,  h=n_per_day)
fc_drift <- forecast(m_drift, h=n_per_day)

accuracy(fc_base,  test)   # MAE, RMSE, MASE…
accuracy(fc_drift, test)

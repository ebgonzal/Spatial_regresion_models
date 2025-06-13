
names(esc_final)
str(esc_final)
summary(esc_final)

vars_resp <- c("I_02", "I_36", "I_710", "I_1114", "I_06", "I_714", "I_314")

# Histograma para cada índice de caries
par(mfrow = c(3, 3))
for (var in vars_resp) {
  hist(esc_final[[var]], main = var, xlab = "Índice de caries", col = "skyblue")
}
par(mfrow = c(1, 1))

vars_sociales <- c("pob", "vincNBI", "preciom2", "Dist_CESAC", "Cob_est")
cor(esc_final[, vars_sociales], use = "complete.obs")

library(ggplot2)

# Graficar cada índice contra cada variable social
for (resp in vars_resp) {
  for (soc in vars_sociales) {
    p <- ggplot(esc_final, aes_string(x = soc, y = resp)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      labs(title = paste("Relación entre", resp, "y", soc)) +
      theme_minimal()
    print(p)
  }
}

pdf("relaciones_indices_variables_sociales.pdf", width = 8, height = 6)

for (resp in vars_resp) {
  for (soc in vars_sociales) {
    p <- ggplot(esc_final, aes_string(x = soc, y = resp)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      labs(title = paste("Relación entre", resp, "y", soc),
           x = soc,
           y = resp) +
      theme_minimal()
    print(p)
  }
}

dev.off()


##############seleccion de mejor modelo GLM##########################
#####################################################################

# Definir función stepwise_forward (asegúrate de tenerla cargada)
stepwise_forward <- function(resp, data, criterio = "AIC", p_threshold = 0.05) {
  todas_vars <- c("Cob_est", "pob", "vincNBI", "preciom2", "Dist_CESAC")
  seleccionadas <- c()
  restantes <- todas_vars
  mejor_criterio <- Inf
  mejor_modelo <- NULL
  historial <- data.frame()
  
  while (length(restantes) > 0) {
    candidatos <- list()
    criterios <- c()
    p_vals <- c()
    
    for (var in restantes) {
      actual <- c(seleccionadas, var)
      f <- as.formula(paste(resp, "~", paste(actual, collapse = " + ")))
      mod <- glm(f, data = data, family = gaussian())
      
      candidatos[[var]] <- mod
      criterios <- c(criterios, ifelse(criterio == "AIC", AIC(mod), BIC(mod)))
      p_vals <- c(p_vals, summary(mod)$coefficients[var, "Pr(>|t|)"])
    }
    
    mejor_idx <- which.min(criterios)
    mejor_var <- restantes[mejor_idx]
    mejor_pval <- p_vals[mejor_idx]
    mejor_mod <- candidatos[[mejor_var]]
    mejor_valor <- criterios[mejor_idx]
    
    if (mejor_valor < mejor_criterio && mejor_pval < p_threshold) {
      seleccionadas <- c(seleccionadas, mejor_var)
      restantes <- setdiff(restantes, mejor_var)
      mejor_criterio <- mejor_valor
      mejor_modelo <- mejor_mod
      historial <- rbind(historial, data.frame(Variable = mejor_var,
                                               AIC = AIC(mejor_mod),
                                               BIC = BIC(mejor_mod),
                                               p_value = mejor_pval))
    } else {
      break
    }
  }
  
  return(list(variables = seleccionadas,
              modelo = mejor_modelo,
              resumen = historial))
}

# Vector con los índices de caries
vars_resp <- c("I_02", "I_36", "I_710", "I_1114", "I_06", "I_714", "I_314")

# Ejecutar stepwise para cada variable de caries
resultados_stepwise <- lapply(vars_resp, function(resp) {
  resultado <- stepwise_forward(resp, esc_final)
  list(
    indice = resp,
    seleccionadas = resultado$variables,
    resumen = resultado$resumen,
    modelo = resultado$modelo
  )
})
names(resultados_stepwise) <- vars_resp

# Mostrar resumen de variables seleccionadas por índice
for (res in resultados_stepwise) {
  cat("\n==============================\n")
  cat("Índice:", res$indice, "\n")
  cat("Variables seleccionadas:", paste(res$seleccionadas, collapse = ", "), "\n")
  print(res$resumen)
}



# Crear una tabla resumen combinada para todos los modelos
tabla_resumen <- do.call(rbind, lapply(resultados_stepwise, function(res) {
  if (nrow(res$resumen) == 0) {
    return(data.frame(
      Indice = res$indice,
      Paso = NA,
      Variable = NA,
      AIC = NA,
      BIC = NA,
      p_value = NA
    ))
  } else {
    res$resumen$Indice <- res$indice
    res$resumen$Paso <- seq_len(nrow(res$resumen))
    res$resumen[, c("Indice", "Paso", "Variable", "AIC", "BIC", "p_value")]
  }
}))

#write.csv(tabla_resumen, "stepwise_resultados.csv", row.names = FALSE)

########################################################################
########################## GWR #########################################


# Guardar modelos seleccionados con nombres como modelo_I_02, modelo_I_36, etc.
for (res in resultados_stepwise) {
  nombre_modelo <- paste0("modelo_", res$indice)
  assign(nombre_modelo, res$modelo, envir = .GlobalEnv)
}

###################################################################
########################## MODELOS GWR ############################

# Convertir a sf y luego a Spatial para GWR y Moran
escuelas_sf <- st_as_sf(esc_final, coords = c("POINT_X", "POINT_Y"), crs = 4326)
escuelas_sp <- as(escuelas_sf, "Spatial")

# Vector de nombres de modelos y sus índices
indices <- c("I_02", "I_36", "I_710", "I_1114", "I_06", "I_714", "I_314")

# Lista para guardar resultados
resultados_moran <- list()
modelos_gwr <- list()

for (indice in indices) {
  
  cat("\n==============================\n")
  cat("📌 Índice:", indice, "\n")
  
  # Acceder al modelo por nombre
  modelo <- get(paste0("modelo_", indice))
  
  # Obtener residuos
  res <- residuals(modelo)
  
  # Crear pesos espaciales
  coords <- coordinates(escuelas_sp)
  vecinos <- knearneigh(coords, k = 5)
  listw <- nb2listw(knn2nb(vecinos))
  
  # Test de Moran
  moran <- moran.test(res, listw)
  resultados_moran[[indice]] <- moran
  
  cat("Moran’s I:", moran$estimate["Moran I statistic"], "\n")
  cat("p-valor:", moran$p.value, "\n")
  
  # Evaluar si hay autocorrelación significativa
  if (moran$p.value < 0.05) {
    cat("⚠️ Autocorrelación detectada. Corriendo GWR...\n")
    
    # Armar fórmula del modelo
    formula_modelo <- formula(modelo)
    
    # Calcular el ancho de banda
    bw <- bw.gwr(formula_modelo,
                 data = escuelas_sp,
                 approach = "AICc",
                 kernel = "bisquare",
                 adaptive = TRUE)
    
    # Ajustar modelo GWR
    gwr_model <- gwr.basic(formula_modelo,
                           data = escuelas_sp,
                           bw = bw,
                           kernel = "bisquare",
                           adaptive = TRUE)
    
    modelos_gwr[[indice]] <- gwr_model
    
    # Exportar coeficientes locales
    gwr_sf <- st_as_sf(gwr_model$SDF)
    st_write(gwr_sf, paste0("gwr_coef_", indice, ".shp"), delete_layer = TRUE)
    
    # Graficar el primer coeficiente local
    primer_var <- names(coef(modelo))[-1][1]
    
    ggplot() +
      geom_sf(data = gwr_sf, aes_string(color = primer_var), inherit.aes = FALSE) +
      scale_color_viridis_c() +
      labs(title = paste("Coeficiente local de", primer_var, "para", indice)) +
      theme_minimal()
    
  } else {
    cat("✅ Sin autocorrelación espacial. GWR no necesario.\n")
  }
}



summary(modelo_I_314)$adj.r.squared
gwr_modelo <- modelos_gwr[["I_314"]]
gwr_modelo$GW.diagnostic$adjR2
gwr_modelo$GW.diagnostic$AICc


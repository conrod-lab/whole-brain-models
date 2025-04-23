library(shiny)
library(cmdstanr)
library(dplyr)
library(ggplot2)
library(posterior)
library(purrr)
library(lme4)
library(lmerTest)
library(cowplot)
library(scales)
library(reshape2)

mod <- cmdstan_model("sim_cortical_thickness.stan")

ui <- fluidPage(
  titlePanel("Cortical Thickness Simulation & ROI Modeling"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n_subj", "Number of Subjects", 100),
      numericInput("n_roi", "Number of ROIs", 200),
      numericInput("n_visit", "Number of Visits", 3),
      
      sliderInput("rho_intra", "Within-network ROI correlation", min = 0, max = 1, value = 0.6),
      sliderInput("rho_inter", "Between-network ROI correlation", min = 0, max = 1, value = 0.2),
      sliderInput("rho_visit", "Temporal correlation across visits", min = 0, max = 1, value = 0.3),
      
      sliderInput("mu_i_sd", "Subject intercept SD", min = 0, max = 1, value = 0.2),
      sliderInput("tau_i_sd", "Subject slope SD", min = 0, max = 1, value = 0.05),
      sliderInput("v_r_sd", "ROI intercept SD", min = 0, max = 1, value = 0.08),
      sliderInput("subj_roi_sd", "Subject × ROI jitter SD", min = 0, max = 1, value = 0.05),
      sliderInput("sigma_eps", "Residual noise SD", min = 0, max = 1, value = 0.4),
      
      sliderInput("gamma_time", "Global Ageing Effect", min = -0.2, max = 0.2, value = -0.1),
      sliderInput("gamma_drug_int", "Drug Main-Effect on Intercept", min = -0.5, max = 0.5, value = 0.2),
      sliderInput("gamma_global", "Drug × Time Interaction", min = -0.5, max = 0.5, value = -0.05),
      
      actionButton("go", "Simulate & Analyze")
    ),
    mainPanel(
      plotOutput("resultsPlot"),
      plotOutput("effectPlot"),
      verbatimTextOutput("console")  
    )
  )
)

server <- function(input, output) {
  
  observeEvent(input$go, {
    
    stan_data <- list(
      n_subj         = input$n_subj,
      n_roi          = input$n_roi,
      n_visit        = input$n_visit,
      n_net          = 5,
      roi_net        = sample(1:5, input$n_roi, replace=TRUE),
      rho_intra      = input$rho_intra,
      rho_inter      = input$rho_inter,
      rho_visit      = input$rho_visit,
      mu_i_sd        = input$mu_i_sd,
      tau_i_sd       = input$tau_i_sd,
      v_r_sd         = input$v_r_sd,
      subj_roi_sd    = input$subj_roi_sd,
      sigma_eps      = input$sigma_eps,
      gamma_time     = input$gamma_time,
      gamma_drug_int = input$gamma_drug_int,
      gamma_global   = input$gamma_global
    )
    withProgress(message = 'Running simulation...', value = 0.5, {
    fit <- mod$sample(data = stan_data, chains = 1, iter_sampling = 1, iter_warmup = 0, fixed_param = TRUE, refresh = 0)
    })
    draws <- fit$draws(format="draws_array")
    
    drug_vec <- as.integer(draws[1,1, grep("^drug\\[", dimnames(draws)$variable)])
    Y_vec <- as.numeric(draws[1,1, grep("^Y\\[", dimnames(draws)$variable)])
    beta_vec <- as.numeric(draws[1,1, grep("^beta_vr\\[", dimnames(draws)$variable)])
    
    Y_array <- array(Y_vec, dim = c(stan_data$n_subj, stan_data$n_visit, stan_data$n_roi))
    beta_mat <- array(beta_vec, dim = c(stan_data$n_visit, stan_data$n_roi))
    
    df <- melt(Y_array, varnames = c("subj", "visit", "roi"), value.name = "thk") %>%
      transmute(
        subj = as.integer(subj),
        visit = as.integer(visit),
        roi = as.integer(roi),
        thk = thk,
        drug = drug_vec[subj]
      ) %>%
      mutate(visit_c = visit - mean(visit), drug_num = drug)
    
    true_df <- tibble(ROI = 1:stan_data$n_roi, true_slope = stan_data$gamma_global + beta_mat[stan_data$n_visit, ])
    
    fit_hier <- lmer(thk ~ visit_c * drug_num + (1 | subj) + (1 + visit_c:drug_num || roi), data = df, REML = FALSE)
    beta_fixed <- fixef(fit_hier)["visit_c:drug_num"]
    coef_tbl <- summary(fit_hier)$coefficients
    p_global <- coef_tbl["visit_c:drug_num", "Pr(>|t|)"]
    sig_pooled <- (p_global < 0.05)
    pooled_df <- ranef(fit_hier)$roi %>%
      as_tibble(rownames="ROI") %>%
      transmute(ROI = as.integer(ROI), slope_pooled = `visit_c:drug_num` + beta_fixed)
    
    unpooled <- map_dfr(1:stan_data$n_roi, function(r) {
      d <- filter(df, roi == r)
      m <- tryCatch(lmer(thk ~ visit_c * drug_num + (1|subj), data = d), error = function(e) NULL)
      if (is.null(m)) return(tibble(ROI=r, slope_unp=NA, p_unp=NA))
      co <- summary(m)$coefficients
      if (!"visit_c:drug_num" %in% rownames(co)) return(tibble(ROI=r, slope_unp=NA, p_unp=NA))
      tibble(ROI = r, slope_unp = co["visit_c:drug_num", "Estimate"], p_unp = co["visit_c:drug_num", "Pr(>|t|)"])
    }) %>% mutate(fdr = p.adjust(p_unp, method="fdr"), sig_fdr = fdr < 0.05)
    
    compare_df <- pooled_df %>%
      left_join(unpooled, by="ROI") %>%
      left_join(true_df, by="ROI")
    
    det_raw <- mean(compare_df$p_unp < 0.05, na.rm=TRUE)
    det_fdr <- mean(compare_df$sig_fdr, na.rm=TRUE)
    rmse <- function(x,y) sqrt(mean((x - y)^2, na.rm=TRUE))
    rmse_unp <- rmse(compare_df$slope_unp, compare_df$true_slope)
    rmse_pld <- rmse(compare_df$slope_pooled, compare_df$true_slope)
    
    p1 <- ggplot(compare_df, aes(x=true_slope, y=slope_unp, color=sig_fdr)) +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="grey") +
      geom_point() +
      scale_color_manual(values=c("grey80","forestgreen")) +
      labs(title="Unpooled vs True ROI Slopes", color="sig") +
      theme_minimal()
    
    det_df <- tibble(Method = c("Pooled", "Unpooled FDR", "Unpooled p<.05"),
                     Detected = c(1, det_fdr, det_raw))
    
    p2 <- ggplot(det_df, aes(Method, Detected, fill=Method)) +
      geom_col() +
      scale_y_continuous(labels=percent_format(accuracy=1), limits=c(0,1)) +
      labs(title="% Effects Detected", y="% tests significant") +
      theme_minimal()
    
    rmse_df <- tibble(Method = c("Pooled", "Unpooled"),
                      RMSE = c(rmse_pld, rmse_unp))
    
    p3 <- ggplot(rmse_df, aes(Method, RMSE, fill=Method)) +
      geom_col() +
      labs(title="RMSE of slope estimates") +
      theme_minimal()
    
    eff_df <- tibble(Method = c("Pooled", "Unpooled Mean", "True"),
                     Estimate = c(beta_fixed, mean(compare_df$slope_unp, na.rm=TRUE), stan_data$gamma_global))
    
    p4 <- ggplot(eff_df, aes(Method, Estimate, fill=Method)) +
      geom_col(data = filter(eff_df, Method != "True"), width = 0.6) +
      geom_point(data = filter(eff_df, Method == "True"),
                 aes(x = Method, y = Estimate),
                 color = "black", size = 4, shape = 18) +
      labs(title = "Recovery of Global Drug×Time Effect", y = "Estimate of γ") +
      theme_minimal()
    
    output$resultsPlot <- renderPlot({
      plot_grid(p1, p2, p3, ncol = 1, rel_heights = c(2, 1, 1))
    })
    
    output$effectPlot <- renderPlot({
      p4
    })
    
    output$console <- renderPrint({
      cat("Simulation Diagnostics:\n")
      cat(sprintf("Pooled model global test p = %.3f → %s\n",
                  p_global,
                  if (sig_pooled) "SIGNIFICANT" else "not significant"))
      
      cat(sprintf("\nUnpooled detection:\n • raw p<.05: %.1f%%\n • FDR<.05: %.1f%%\n",
                  det_raw * 100, det_fdr * 100))
      
      cat(sprintf("\nRMSE:\n • Unpooled: %.3f\n • Pooled:   %.3f\n",
                  rmse_unp, rmse_pld))
      
      cat(sprintf("\nGlobal effect recovery:\n"))
      cat(sprintf("  True γ:     %.3f\n", stan_data$gamma_global))
      cat(sprintf("  Pooled:     %.3f  (Bias %+ .3f)\n",
                  beta_fixed, beta_fixed - stan_data$gamma_global))
      cat(sprintf("  Unpooled:   %.3f  (Bias %+ .3f)\n\n",
                  mean(compare_df$slope_unp, na.rm = TRUE),
                  mean(compare_df$slope_unp, na.rm = TRUE) - stan_data$gamma_global))
      
      cat("=== Summary of Pooled Mixed Model ===\n")
      print(summary(fit_hier))
      
      # Optional: Show one unpooled model summary
      example_model <- tryCatch(
        lmer(thk ~ visit_c * drug_num + (1|subj), data = filter(df, roi == 1)),
        error = function(e) NULL
      )
      if (!is.null(example_model)) {
        cat("\n=== Summary of Unpooled Model for ROI 1 ===\n")
        print(summary(example_model))
      } else {
        cat("\n(Unpooled model for ROI 1 failed to converge or fit.)\n")
      }
    })
    
    
  })
}

shinyApp(ui = ui, server = server)

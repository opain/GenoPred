#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(cowplot)

# Create function
which_quant <- function(PRS_z_score=0, n_quantile=20){
    E_PRS = 0
    SD_PRS = sqrt(1)

    by_quant<-1/(n_quantile)
    PRS_quantile_bounds <- qnorm(p=seq(0, 1, by=by_quant))
    lower_PRS_vec <- PRS_quantile_bounds[1:n_quantile]
    upper_PRS_vec <- PRS_quantile_bounds[2:(n_quantile+1)]
    
    out<-data.frame(q=1:n_quantile,
                    q_min=lower_PRS_vec,
                    q_max=upper_PRS_vec)

    quant<-out$q[PRS_z_score > out$q_min & PRS_z_score <= out$q_max]
    return(quant)
}
    
mean_sd_quant.f <- function(PRS_R2=0.641, Outcome_mean=1, Outcome_sd=1, n_quantile=20, quant=NA){
    ### PRS quantiles with a continuous phenotype (Y)
    library(tmvtnorm)
    ###
    E_PRS = 0
    SD_PRS = sqrt(1)
    E_phenotype = Outcome_mean
    SD_phenotype = Outcome_sd 

    by_quant<-1/(n_quantile)
    PRS_quantile_bounds <- qnorm(p=seq(0, 1, by=by_quant), mean= E_PRS, sd= SD_PRS)
    lower_PRS_vec <- PRS_quantile_bounds[1:n_quantile]
    upper_PRS_vec <- PRS_quantile_bounds[2:(n_quantile+1)]
    
    mean_vec <- c(E_phenotype, E_PRS)
    sigma_mat <- matrix(sqrt(PRS_R2)*SD_phenotype*SD_PRS, nrow=2, ncol=2)
    sigma_mat[1,1] <- SD_phenotype^2
    sigma_mat[2,2] <- SD_PRS^2
    
    ### mean of phenotype within the truncated PRS distribution
    out_mean_Y <- rep(0, 20)
    ### SD of phenotype within the truncated PRS distribution
    out_SD_Y <- rep(0, 20)
    ### cov of Y and PRS given truncation on PRS
    out_cov_Y_PRS <- rep(0, 20)
    ### SD of PRS given truncation on PRS
    out_SD_PRS <- rep(0, 20)
    ### mean PRS given truncation on PRS
    out_mean_PRS <- rep(0, 20)
    
    if(!is.na(quant)){
        i<-quant
        
        distribution_i <- mtmvnorm(mean = mean_vec,
                                   sigma = sigma_mat,
                                   lower = c(-Inf, lower_PRS_vec[i]),
                                   upper = c(Inf, upper_PRS_vec[i]),
                                   doComputeVariance=TRUE,
                                   pmvnorm.algorithm=GenzBretz())
        out_mean_Y[i] <- distribution_i$tmean[1]
        out_mean_PRS[i] <- distribution_i$tmean[2]
        out_SD_Y[i] <- sqrt(distribution_i$tvar[1,1])
        out_SD_PRS[i] <- sqrt(distribution_i$tvar[2,2])
        out_cov_Y_PRS[i] <- distribution_i$tvar[1,2]
        
        out<-data.frame(q=quant,
                        q_min=lower_PRS_vec[quant],
                        q_max=upper_PRS_vec[quant],
                        x_mean=out_mean_Y[quant],
                        x_sd=out_SD_Y[quant])
        
    } else {
        for(i in 1:n_quantile){
            distribution_i <- mtmvnorm(mean = mean_vec,
                                       sigma = sigma_mat,
                                       lower = c(-Inf, lower_PRS_vec[i]),
                                       upper = c(Inf, upper_PRS_vec[i]),
                                       doComputeVariance=TRUE,
                                       pmvnorm.algorithm=GenzBretz())
            out_mean_Y[i] <- distribution_i$tmean[1]
            out_mean_PRS[i] <- distribution_i$tmean[2]
            out_SD_Y[i] <- sqrt(distribution_i$tvar[1,1])
            out_SD_PRS[i] <- sqrt(distribution_i$tvar[2,2])
            out_cov_Y_PRS[i] <- distribution_i$tvar[1,2]
        }
        
        out<-data.frame(q=1:n_quantile,
                        q_min=lower_PRS_vec,
                        q_max=upper_PRS_vec,
                        x_mean=out_mean_Y,
                        x_sd=out_SD_Y)
    }
    
    return(out)
    
    out_mean_Y
    out_SD_Y
    
    out_mean_PRS
    out_SD_PRS
    out_cov_Y_PRS
}

# Define UI for application
ui <- fluidPage(
    
    # Application title
    titlePanel("Polygenic Score"),
    
    # Sidebar with a slider input for PRS z-score
    sidebarLayout(
        sidebarPanel(
            sliderInput("z_score",
                        "Polygenic Z-score:",
                        min = -3,
                        max = 3,
                        value = 0,
                        step=0.01),

            sliderInput("r2",
                        "Variance explained by polygenic score (%):",
                        min = 0.01,
                        max = 99.9,
                        value = 30,
                        step=0.1),
            
            numericInput("pop_mean", label = "Population mean of trait:", value = 0),

            numericInput("pop_sd", label = "Population standard deviation of trait:", value = 1),
            
            sliderInput("ci",
                        "Prediction confidence interval (%):",
                        min = 50,
                        max = 100,
                        value = 95,
                        step=0.1)
            
        ),
        
        # Show plot
        mainPanel(
            h3("Relative Risk"),
            plotOutput("Plot_rel", height = "325px", width = "600px"),
            h3("Absolute Risk"),
            plotOutput("Plot_abs", height = "350px", width = "600px")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$Plot_rel <- renderPlot({
        PRS_z_score<-input$z_score

        prs_dist<- ggplot(data = data.frame(x = c(-4, 4)), aes(x=x)) +
            stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1)) +
            stat_function(fun = dnorm, args = list(mean = 0, sd = 1), xlim = c(PRS_z_score, -4),
                          geom = "area", fill = "#CC66FF", alpha = .4) +
            stat_function(fun = dnorm, args = list(mean = 0, sd = 1), xlim = c(PRS_z_score, 4),
                          geom = "area", fill = "#FF6633", alpha = .4) +
            geom_vline(xintercept=PRS_z_score, linetype='dashed') +
            geom_text(label=paste0(round(pnorm(PRS_z_score)*100,1),"% have lower \npolygenic scores"), mapping=aes(x=PRS_z_score-0.1, y=0.5), colour='#CC66FF', hjust='right', vjust=0.8, size=5) +
            geom_text(label=paste0(round(100-(pnorm(PRS_z_score)*100),1),"% have higher \npolygenic scores"), mapping=aes(x=PRS_z_score+0.1, y=0.5), colour='#FF6633', hjust='left', vjust=0.8, size=5) +
            scale_y_continuous(breaks = NULL) +
            theme_half_open() +
            xlim(-5,5) +
            labs(y='Number of people', x='Polygenic Score', title='Distribution of polygenic scores') +
            theme(plot.title = element_text(hjust = 0.5))
        
        plot_grid(prs_dist, ncol = 1)
        
    })
    
    output$Plot_abs <- renderPlot({
        # Define parameters
        PRS_z_score<-input$z_score
        PRS_R2=input$r2/100
        Outcome_mean=input$pop_mean
        Outcome_sd=input$pop_sd
        conf_int=input$ci/100
        n_quant<-1000
        
        # Run analysis
        quant<-which_quant(PRS_z_score = PRS_z_score, n_quantile = n_quant)
        
        risk_quantiles<-mean_sd_quant.f(PRS_R2=PRS_R2, Outcome_mean=Outcome_mean, Outcome_sd=Outcome_sd, n_quantile=n_quant, quant=quant)
        indiv_result_all<-risk_quantiles[PRS_z_score > risk_quantiles$q_min & PRS_z_score <= risk_quantiles$q_max,]
        
        ref<-dnorm(seq(-4.5,4.5,0.01), 0, 1)
        ref_plot<-data.frame(x=(seq(-4.5,4.5,0.01)*Outcome_sd)+Outcome_mean,
                             y=ref,
                             Group='General public')
        
        indiv<-dnorm(seq(-4.5,4.5,0.01), 0, 1)
        indiv_plot<-data.frame(x=(seq(-4.5,4.5,0.01)*indiv_result_all$x_sd)+indiv_result_all$x_mean,
                               y=indiv,
                               Group='People like you')
        
        plot_dat<-rbind(ref_plot, indiv_plot)
        plot_dat$Group<-factor(plot_dat$Group, levels=c('General public','People like you'))
        
        # Calculate 95CI for target individual
        lowCI<-qnorm((1-conf_int)/2,indiv_result_all$x_mean,indiv_result_all$x_sd)
        highCI<-qnorm(1-((1-conf_int)/2),indiv_result_all$x_mean,indiv_result_all$x_sd)
        
        abs_dist<-ggplot(plot_dat, aes(x=x, y=y, fill=Group)) +
            geom_area(alpha=0.4, colour='black') +
            scale_fill_manual(values=c("#84CA72","#0066CC")) +
            labs(y='Number of people', x='Trait', title='You compared to general population', fill=NULL) +
            geom_segment(aes(x = indiv_result_all$x_mean , y = 0, xend = indiv_result_all$x_mean, yend = 0.4), color="black") +
            geom_segment(aes(x = lowCI , y = 0, xend = lowCI, yend = 0.38), color="#0066CC", linetype="dashed") +
            geom_segment(aes(x = highCI , y = 0, xend = highCI, yend = 0.38), color="#0066CC", linetype="dashed") +
            geom_text(label=paste0('Estimate = ',round(indiv_result_all$x_mean,2)), mapping=aes(x=Outcome_mean+Outcome_sd, y=0.5), colour='black', hjust='left', vjust=0.8, size=5, check_overlap = TRUE) +
            geom_text(label=paste0(conf_int*100,'% CI = ',round(lowCI,2),' – ',round(highCI,2)), mapping=aes(x=Outcome_mean+Outcome_sd, y=0.5), colour='#0066CC', hjust='left', vjust=2.5, size=5, check_overlap = TRUE) +
            scale_y_continuous(breaks = NULL) +
            theme_half_open() +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.position=c(0.01,0.95), legend.box = "horizontal")
        plot_grid(abs_dist, ncol = 1)
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

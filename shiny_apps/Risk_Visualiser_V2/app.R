#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

start_time <- Sys.time()

library(shiny)
library(ggplot2)
library(reshape2)
library(cowplot)
library(scales)
library(ggchicklet)

# Create effect size conversion functions

r2_r2l <- function(k, r2l, p) {
    #Lee SH, Goddard ME, Wray NR, Visscher PM. (2012)
    x= qnorm(1-k)
    z= dnorm(x)
    i=z/k
    C= k*(1-k)*k*(1-k)/(z^2*p*(1-p))
    theta= i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
    r = sqrt(r2l)/sqrt(C-(C*(r2l)*theta))
    r2<-r^2
    r2
}

d_r2<-function(r2,p){
    r<-sqrt(r2)
    
    n_case<-p
    n_con<-1-p
    
    a<-(n_case+n_con)^2/(n_case*n_con)
    
    d<-sqrt(a)*r/sqrt(1-r^2)
    d
}

# Create function for z-score to absolute risk
ccprobs.f <- function(d=0.641, prev=0.7463, n_quantile=20){
    mu_case <- d
    mu_control <- 0
    
    varPRS <- prev*(1+(d^2) - (d*prev)^2) + (1-prev)*(1 - (d*prev)^2)
    E_PRS <- d*prev
    
    by_quant<-1/n_quantile
    p_quant <- seq(by_quant, 1-by_quant, by=by_quant)
    quant_vals_PRS <- rep(0, length(p_quant))
    quant_f_solve <- function(x, prev, d, pq){prev*pnorm(x-d) + (1-prev)*pnorm(x) - pq}
    for(i in 1:length(p_quant)){
        quant_vals_PRS[i] <- unlist(uniroot(quant_f_solve, prev=prev, d=d, pq= p_quant[i], interval=c(-2.5, 2.5), extendInt = "yes", tol=6e-12)$root)
    }
    
    ul_qv_PRS <- matrix(0, ncol=2, nrow=n_quantile)
    ul_qv_PRS[1,1] <- -Inf
    ul_qv_PRS[2:n_quantile,1] <- quant_vals_PRS
    ul_qv_PRS[1:(n_quantile-1),2] <- quant_vals_PRS
    ul_qv_PRS[n_quantile,2] <- Inf
    
    ul_qv_PRS<-cbind(ul_qv_PRS, (ul_qv_PRS[,1:2]-E_PRS)/sqrt(varPRS))
    
    prob_quantile_case <- pnorm(ul_qv_PRS[,2], mean = mu_case) - pnorm(ul_qv_PRS[,1], mean = mu_case)
    prob_quantile_control <- pnorm(ul_qv_PRS[,2], mean = mu_control) - pnorm(ul_qv_PRS[,1], mean = mu_control)
    p_case_quantile <- (prob_quantile_case*prev)/by_quant
    p_cont_quantile <- (prob_quantile_control*(1-prev))/by_quant
    
    OR <- p_case_quantile/p_cont_quantile
    OR <- OR/OR[1]
    out <- cbind(ul_qv_PRS[,3:4],p_cont_quantile, p_case_quantile, OR)
    row.names(out) <- 1:n_quantile
    colnames(out) <- c("q_min", "q_max","p_control", "p_case", "OR")
    
    data.frame(out)
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
            
            br(),
            
            sliderInput("prev",
                        "Population Prevalence (%):",
                        min = 0.1,
                        max = 99.9,
                        value = 50,
                        step=0.1),

            br(),
            
            sliderInput("samp",
                        "Sampling Fraction (%):",
                        min = 0.1,
                        max = 99.9,
                        value = 50,
                        step=0.1),
            
            br(),
            
            radioButtons("type", "Polygenic Score Effect Size Type:",
                         c("Cohen's D" = "d",
                           "Area Under-the-ROC Curve (AUC)" = "auc",
                           "Odds Ratio (1SD)" = "or",
                           "Observed R-squared" = "r2_obs",
                           "Liability R-squared" = "r2_liab")),
            
            numericInput("val", "Polygenic Score Effect Size Value", 0.7, step = 0.001),
            
            br(),
        ),
        
        # Show plot
        mainPanel(
            plotOutput("Plot", height = "1400px", width = "600px")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$Plot <- renderPlot({
        k<-input$prev/100
        p<-input$samp/100
        
        # Convert effect size to cohen's d
        if(input$type == 'd'){
            d<-input$val
        }
        if(input$type == 'auc'){
            auc<-input$val
            d<-sqrt(2)*qnorm(auc)
        }
        if(input$type == 'or'){
            or<-input$val
            f<-function(d,OR,p){OR - (exp(d*sqrt(1 + d^2*p*(1-p))))}
            d<-uniroot(f, p=p, OR=or, interval=c(-1, 1), extendInt = "yes", tol=6e-12)$root
        }
        if(input$type == 'r2_obs'){
            r2_obs<-input$val
            d<-d_r2(r2=r2_obs, p=p)
        }
        if(input$type == 'r2_liab'){
            r2_liab<-input$val
            r2<-r2_r2l(k=k,r2l=r2_liab,p=p)
            d<-d_r2(r2=r2, p=p)
        }
        
        # Define parameters
        trait<-'trait'
        PRS_z_score<-input$z_score
        prev<-input$prev/100

        risk_quantiles<-ccprobs.f(d=d, prev=prev, n_quantile=1000)
        
        indiv_result_all<-risk_quantiles[PRS_z_score > risk_quantiles$q_min & PRS_z_score < risk_quantiles$q_max,]
        indiv_result<-indiv_result_all[, c('p_case', 'p_control')]
        
        indiv_result<-melt(indiv_result)
        indiv_result$variable<-c('Case','Control')
        indiv_result$variable<-factor(indiv_result$variable,levels = c('Control','Case'))
        
        # Create a grid of dots, some red and some blue
        n_case_1<-round(indiv_result$value*100,1)[1]
        n_control_1<-round(indiv_result$value*100,1)[2]
        n_case<-round(indiv_result$value*100)[1]
        n_control<-round(indiv_result$value*100)[2]
        
        plot_dat<-(matrix(c(rep('Case',n_case),rep('Control',n_control)),nrow=10, ncol=10))
        plot_dat<-data.frame(melt(plot_dat))
        plot_dat$value<-factor(plot_dat$value, levels = c('Control','Case'))
        
        # Create data for general population example
        pop_case_1<-round(input$prev,1)
        pop_case<-round(input$prev)
        pop_control_1<-round(100-input$prev,1)
        pop_control<-round(100-input$prev)
        pop_result<-data.frame(variable=c('Case','Control'),
                               value=c(prev,1-prev))
        
        pop_result$variable<-factor(pop_result$variable, levels = c('Control','Case'))
        
        plot_dat_pop<-(matrix(c(rep('Case',pop_case),rep('Control',100-pop_case)),nrow=10, ncol=10))
        plot_dat_pop<-data.frame(melt(plot_dat_pop))
        plot_dat_pop$value<-factor(plot_dat_pop$value, levels = c('Control','Case'))
        
        dnorm_new<-function(x, mean=0, sd=1, log=F, height=1){
            out<-dnorm(x=x,mean=mean,sd=sd,log=log)
            out*height
        }
        
        dnorm_2_new<-function(x, mean_1=0, mean_2=0, sd_1=1, sd_2=1, log=F, p_2=0.5){
            out_1<-dnorm_new(x, mean=mean_1, sd=sd_1, log=log, height=1-p_2)
            out_2<-dnorm_new(x, mean=mean_2, sd=sd_2, log=log, height=p_2)
            out_1+out_2
        }

        varPRS <- prev*(1+(d^2) - (d*prev)^2) + (1-prev)*(1 - (d*prev)^2)
        sdPRS<-sqrt(varPRS)
        E_PRS <- d*prev
        
        prs_dist<- ggplot(data = data.frame(x = c(-4, 4)), aes(x=x)) +
            stat_function(fun = dnorm_2_new, n = 101, args = list(mean_1 = 0-(E_PRS), sd_1 = 1, mean_2 = d-(E_PRS), sd_2 = 1, p_2=prev)) +
            ylab("") +
            stat_function(fun = dnorm_2_new, args = list(mean_1 = 0-(E_PRS), sd_1 = 1, mean_2 = d-(E_PRS), sd_2 = 1, p_2=prev), xlim = c(PRS_z_score, -5),
                          geom = "area", fill = "#84CA72", alpha = .4) +
            stat_function(fun = dnorm_2_new, args = list(mean_1 = 0-(E_PRS), sd_1 = 1, mean_2 = d-(E_PRS), sd_2 = 1, p_2=prev), xlim = c(PRS_z_score, 5),
                          geom = "area", fill = "#0066CC", alpha = .4) +
            geom_vline(xintercept=PRS_z_score, linetype='dashed') +
            geom_text(label=paste0(round(pnorm(PRS_z_score)*100,1),"% have lower \npolygenic scores"), mapping=aes(x=PRS_z_score-0.1, y=0.5), colour='#84CA72', hjust='right', vjust=0.8, size=5) +
            geom_text(label=paste0(round(100-(pnorm(PRS_z_score)*100),1),"% have higher \npolygenic scores"), mapping=aes(x=PRS_z_score+0.1, y=0.5), colour='#0066CC', hjust='left', vjust=0.8, size=5) +
            scale_y_continuous(breaks = NULL) +
            theme_half_open() +
            xlim(-5,5) +
            labs(y='Number of people', x='Polygenic Score', title='Distribution of polygenic scores') +
            theme(plot.title = element_text(hjust = 0.5))
        
        prs_dist_2<- ggplot(data = data.frame(x = c(-4, 4)), aes(x=x)) +
            stat_function(fun = dnorm_2_new, n = 101, args = list(mean_1 = 0-(E_PRS), sd_1 = 1, mean_2 = d-(E_PRS), sd_2 = 1, p_2=prev)) +
            stat_function(fun = dnorm_new, n = 101, args = list(mean = 0-(E_PRS), sd = 1, height=1-prev)) +
            stat_function(fun = dnorm_new, n = 101, args = list(mean = d-(E_PRS), sd = 1, height=prev)) +
            stat_function(fun = dnorm_new, args = list(mean = 0-(E_PRS), sd = 1, height=1-prev), geom = "area", fill = "#84CA72", alpha = .4) +
            stat_function(fun = dnorm_new, args = list(mean = d-(E_PRS), sd = 1, height=prev), geom = "area", fill = "#0066CC", alpha = .4) +
            geom_vline(xintercept=PRS_z_score, linetype='dashed') +
            geom_text(label=paste0(n_case_1,'% of people\nlike you have ',trait), mapping=aes(x=PRS_z_score+0.1, y=0.5), colour='#0066CC', hjust='left', vjust=0.8, size=5) +
            scale_y_continuous(breaks = NULL) +
            theme_half_open() +
            xlim(-5,5) +
            labs(y='Number of people', x='Polygenic Score', title='Distribution of polygenic scores') +
            theme(plot.title = element_text(hjust = 0.5))
        
        bar_chart<-ggplot(data.frame(x=1,y=0:1), aes(x=x, y=y)) +
            geom_chicklet(radius = grid::unit(1, 'mm'), data=indiv_result, mapping=aes(x=1, y=value, fill=variable), stat="identity",position='stack') +
            scale_fill_manual(values=c("#84CA72","#0066CC"), drop = F) +
            annotate("text", x=1.5, y=((((n_control_1)/2))+n_case_1)/100, label=paste0(n_control_1,"%\ndo not have ",trait), colour = '#84CA72', hjust=0, size=6) +
            annotate("text", x=1.5, y=((n_case_1/2))/100, label=paste0(n_case_1,'%\nhave ',trait), colour = '#0066CC', hjust=0, size=6) +
            ylim(-0.1,1.05) +
            theme_half_open() +
            labs(title='Of people with your genetics,') +
            xlim(0.25,2.5) +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none")
        
        bar_chart_pop<-ggplot(data.frame(x=1,y=0:1), aes(x=x, y=y)) +
            geom_chicklet(radius = grid::unit(1, 'mm'), data=pop_result, mapping=aes(x=1, y=value, fill=variable), stat="identity",position='stack') +
            scale_fill_manual(values=c("#84CA72","#0066CC"), drop = F) +
            annotate("text", x=1.5, y=(((pop_control_1/2))+pop_case_1)/100, label=paste0(pop_control_1,"%\ndo not have ",trait), colour = '#84CA72', hjust=0, size=6) +
            annotate("text", x=1.5, y=((pop_case_1/2))/100, label=paste0(pop_case_1,'%\nhave ',trait), colour = '#0066CC', hjust=0, size=6) +
            ylim(-0.1,1.05) +
            theme_half_open() +
            labs(title='In the general population,') +
            xlim(0.25,2.5) +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none")
        
        bar_chart_grid<-plot_grid(bar_chart, bar_chart_pop, labels = NULL, nrow = 1)
        
        person_plot<-
            ggplot(plot_dat, aes(x=Var1, y=Var2, colour=value)) +
            geom_point(mapping=aes(x=Var1, y=Var2, colour=value), size=2) +
            geom_point(mapping=aes(x=Var1, y=Var2-0.3, colour=value), size=4.5) +
            scale_colour_manual(values=c("#84CA72","#0066CC"), drop = F) +
            annotate("text", x=12, y=((n_case/2))/10+0.5, label=paste0(n_case,"%\nhave ",trait), colour = '#0066CC', hjust=0, size=6) +
            annotate("text", x=12, y=(((100-n_case)/2)+n_case)/10, label=paste0(100-n_case,"%\ndon't have ",trait), colour = '#84CA72', hjust=0, size=6) +
            xlim(0.5,20) +
            ylim(0,11) +
            theme_half_open() +
            labs(title='Of people with your genetics,') +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none")

        person_plot_pop<-
            ggplot(plot_dat_pop, aes(x=Var1, y=Var2, colour=value)) +
            geom_point(mapping=aes(x=Var1, y=Var2, colour=value), size=2) +
            geom_point(mapping=aes(x=Var1, y=Var2-0.3, colour=value), size=4.5) +
            scale_colour_manual(values=c("#84CA72","#0066CC"), drop = F) +
            annotate("text", x=12, y=((pop_case/2))/10+0.5, label=paste0(pop_case,"%\nhave ",trait), colour = '#0066CC', hjust=0, size=6) +
            annotate("text", x=12, y=(((100-pop_case)/2)+pop_case)/10, label=paste0(100-pop_case,"%\ndon't have ",trait), colour = '#84CA72', hjust=0, size=6) +
            xlim(0.5,20) +
            ylim(0,11) +
            theme_half_open() +
            labs(title='In the general population,') +
            theme(axis.line=element_blank(),axis.text.x=element_blank(),
                  axis.text.y=element_blank(),axis.ticks=element_blank(),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank(),
                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),plot.background=element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none")
        
        person_plot_grid<-plot_grid(person_plot, person_plot_pop, labels = NULL, nrow = 1)
        
        plot_grid(prs_dist, prs_dist_2, bar_chart_grid, person_plot_grid, ncol = 1)
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

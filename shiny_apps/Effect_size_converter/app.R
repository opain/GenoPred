#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)

start_time <- Sys.time()

r2l_r2 <- function(k, r2, p) {
    #Lee SH, Goddard ME, Wray NR, Visscher PM. (2012)
    x= qnorm(1-k)
    z= dnorm(x)
    i=z/k
    C= k*(1-k)*k*(1-k)/(z^2*p*(1-p))
    theta= i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
    h2l_R2 = C*r2 / (1 + C*theta*r2)
    h2l_R2
}

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

r2_d<-function(d,p){
    n_case<-p
    n_con<-1-p
    
    a<-(n_case+n_con)^2/(n_case*n_con)
    
    r=d/sqrt(a+d^2)
    r2<-r^2
    r2
}

library(shiny)

# Define UI for application
ui <- fluidPage(
    fluidRow(
        column(3,
            
            radioButtons("type", "Effect size Metric:",
                         c("Cohen's D" = "d",
                           "Area Under-the-ROC Curve (AUC)" = "auc",
                           "Odds Ratio (1SD)" = "or",
                           "Observed R-squared" = "r2_obs",
                           "Liability R-squared" = "r2_liab"))),
        column(3,
            numericInput("val", "Effect Size Value:", 0.7, step = 0.001),
            numericInput("p", "Sampling Fraction:", 0.5, min = 0, max=1, step = 0.001),
            numericInput("k", "Population Prevalence:", 0.5, min = 0, max=1, step = 0.001))
    ),
        
    # Show plot
    fluidRow(column(6, dataTableOutput('res')))
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    url <- a("https://opain.github.io/GenoPred/", href="https://opain.github.io/GenoPred/")
    output$tab <- renderUI({
        tagList("Description of conversion methodology:", url)
    })
    
    output$res <- renderDataTable({
        if(input$type == 'd'){
            
            d<-input$val
            auc<-pnorm(abs(d)/sqrt(2), 0, 1)
            
            tmp<-data.frame(d=d,
                            auc=auc,
                            or=exp(d*sqrt(1 + d^2*input$p*(1-input$p))),
                            r2_obs=r2_d(d=d,p=input$p),
                            r2_liab=r2l_r2(k=input$k, r2=r2_d(d=d,p=0.5), p=0.5),
                            p=input$p,
                            k=input$k)
            
        }
        if(input$type == 'auc'){
            
            auc<-input$val
            d<-sqrt(2)*qnorm(auc)
            
            tmp<-data.frame(d=d,
                            auc=auc,
                            or=exp(d*sqrt(1 + d^2*input$p*(1-input$p))),
                            r2_obs=r2_d(d=d,p=input$p),
                            r2_liab=r2l_r2(k=input$k, r2=r2_d(d=d,p=0.5), p=0.5),
                            p=input$p,
                            k=input$k)
            
        }
        if(input$type == 'or'){
            
            or<-input$val
            f<-function(d,OR,p){OR - (exp(d*sqrt(1 + d^2*p*(1-p))))}
            d<-uniroot(f, p=input$p, OR=or, interval=c(-1, 1), extendInt = "yes", tol=6e-12)$root
            auc<-pnorm(abs(d)/sqrt(2), 0, 1)
            
            tmp<-data.frame(d=d,
                            auc=auc,
                            or=or,
                            r2_obs=r2_d(d=d,p=input$p),
                            r2_liab=r2l_r2(k=input$k, r2=r2_d(d=d,p=0.5), p=0.5),
                            p=input$p,
                            k=input$k)
            
        }
        if(input$type == 'r2_obs'){
            
            r2_obs<-input$val
            d<-d_r2(r2=r2_obs, p=input$p)
            auc<-pnorm(abs(d)/sqrt(2), 0, 1)
            
            tmp<-data.frame(d=d,
                            auc=auc,
                            or=exp(d*sqrt(1 + d^2*input$p*(1-input$p))),
                            r2_obs=r2_obs,
                            r2_liab=r2l_r2(k=input$k,r2=r2_obs,p=input$p),
                            p=input$p,
                            k=input$k)
            
        }
        if(input$type == 'r2_liab'){
            
            r2_liab<-input$val
            r2<-r2_r2l(k=input$k,r2l=r2_liab,p=input$p)
            d<-d_r2(r2=r2, p=input$p)
            auc<-pnorm(abs(d)/sqrt(2), 0, 1)
            
            tmp<-data.frame(d=d,
                            auc=auc,
                            or=exp(d*sqrt(1 + d^2*input$p*(1-input$p))),
                            r2_obs=r2_d(d=d,p=input$p),
                            r2_liab=r2_liab,
                            p=input$p,
                            k=input$k)
            
        }
        
        names(tmp)<-c("Cohen's <i>d</i>", "AUC", "OR (1SD)","Observed <i>R</i><sup>2</sup>","Liability <i>R</i><sup>2</sup>","Sampling","Prevalence")
        
        datatable(round(tmp,3), escape = FALSE,options = list(dom = 't',ordering=F,columnDefs = list(list(className = 'dt-center', targets = '_all'))), rownames = F)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

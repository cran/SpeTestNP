## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, comment = "#>")

## ----eval = F-----------------------------------------------------------------
#  install.packages("SpeTestNP")
#  

## ----eval = F-----------------------------------------------------------------
#  install.packages("devtools")
#  
#  library("devtools")
#  
#  install_github("HippolyteBoucher/SpeTestNP")
#  

## ----eval = F-----------------------------------------------------------------
#    SpeTest(eq)
#  

## -----------------------------------------------------------------------------

    library(SpeTestNP)
    library(AER)

    ### Loading the data and taking a first look

    data( CPSSW8 )
    
    summary ( CPSSW8 )
    

## -----------------------------------------------------------------------------

    lm_lin <- lm( earnings ~ age + education,
                      data = CPSSW8[1:1000,] )

    summary ( lm_lin )


## -----------------------------------------------------------------------------

    SpeTest( lm_lin , type = "icm" , rejection = "bootstrap" )
    SpeTest( lm_lin , type = "zheng" , rejection = "asymptotics" )


## -----------------------------------------------------------------------------

    lm_quad <- lm( earnings ~ age + I(age^2) + education + I(education^2),
                      data = CPSSW8[1:1000,] )

    summary( lm_quad )
    
    SpeTest( lm_quad , type = "icm" , rejection = "bootstrap" )
    SpeTest( lm_quad , type = "zheng" , rejection = "asymptotics")


## -----------------------------------------------------------------------------

    lm_nlin <- lm( earnings ~ age + I(age^2) + education + I(education^2) 
                       + I(education*age) + I(education^2*age)
                       + I(education*age^2) + I(education^2*age^2),
                       data= CPSSW8[1:1000,] )

    summary( lm_nlin )
    
    SpeTest( lm_nlin , type = "icm" , rejection = "bootstrap" )
    SpeTest( lm_nlin , type = "zheng" , rejection = "asymptotics")


## -----------------------------------------------------------------------------

    SpeTest( lm_nlin, type = "pala", rejection = "asymptotics", nbeta = 40 )
    SpeTest( lm_nlin, type = "pala", rejection = "bootstrap" , nboot = 10 , nbeta = 10 )



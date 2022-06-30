line_plot <- function(op_chars,xvar,yvar,
                      xlab='g',
                      ylab='g',
                      save=FALSE) {
  #' Function for generating plots of operating characteristics
  #'
  #' INPUTS:
  #' - op_chars
  #' - xvar (u,a,n,b)
  #' - yvar (bias,mse,ATE,ci_cov,power)
  #' - xlab
  #' - ylab
  #' - save
  #' - filename
  #'
  #' OUTPUTS:
  #' Nothing, but will produce a plot in the plotting window if using interactively,
  #' and will save a file with name format xvar_yvar_data.pdf if save is set to TRUE
  #' 
  
  # 
  op_chars %>% ggplot(aes(x=eval(as.name(xvar)),
                          y=eval(as.name(yvar)),
                          color=method)) + 
               geom_line() + theme_bw() +
               xlab(xlab) + ylab(ylab)
}

line_plot(op_chars,'u')

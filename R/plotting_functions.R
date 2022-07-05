set_plot_theme <- function() {
  # theme_set(theme_bw() +
  #             theme(plot.title = element_text(hjust = 0, size = 16),
  #                   plot.subtitle = element_text(hjust = 0, size = 12),
  #                   axis.title = element_text(size = 12),
  #                   strip.text = element_text(size = 12),
  #                   legend.position = "bottom"))
  ggplot2::theme_set(
    ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        plot.subtitle = ggplot2::element_text(
          margin = ggplot2::margin(b = 10)
        ),
        strip.text = ggplot2::element_text(
          face = "bold",
          hjust = 0,
          margin = ggplot2::margin(b = 5),
          size = 11
        )
      )
  )
}

line_plot <- function(op_chars,xvar,yvar,
                      xlab,
                      ylab,
                      plt_title,
                      fixed_vars=c('u','a','n','b'),
                      fixed_defaults=c(0.3,1,5000,0),
                      var_labs=c('ME variance:', 'True ATE:', 'n:', 'Binary outcome:'),
                      drop_methods=NULL,
                      extra_options=NULL,
                      save=FALSE) {
  #' Function for generating plots of operating characteristics
  #'
  #' INPUTS:
  #' - op_chars
  #' - xvar (u,a,n,b): A string containing the x variable for the plot
  #' - yvar (bias,mse,ATE,ci_cov,power): A string for y variable
  #' - xlab: x axis title
  #' - ylab: y axis title 
  #' - drop_methods: A vector of strings for methods you would like to not be 
  #'                 included in the plot (e.g. could drop the naive method)
  #' - extra_options: A string containing any additional ggplot options you want
  #'                  - must be written in ggplot syntax
  #' - save: A logical that, when set to TRUE, will save a plot in pdf format
  #'         with the following name convention: <xvar>_<yvar>_<fixedvars/vals>.pdf
  #'         e.g. u_power_n5000_a1_b0.pdf would be the name of a line plot of
  #'         power over varying error variance, with n=5000, the ATE=1 and the 
  #'         outcome being continuous
  #'
  #' OUTPUTS:
  #' Nothing, but will produce a plot in the plotting window if using interactively,
  #' and will save a file if save is set to TRUE
  #' 
  
  # Use the custom plot theme
  set_plot_theme()
  
  # If don't want to include any methods, filter them out
  for (m in drop_methods) {
    op_chars <- op_chars %>% filter(method!=m)
  }
  
  # Prune out the var that's being varied
  drop_idx <- (fixed_vars %in% xvar)   
  fixed_vars <- fixed_vars[! drop_idx] ; fixed_defaults <- fixed_defaults[! drop_idx]
  var_labs <- var_labs[! drop_idx]
  
  
  # Fix values at defaults via filtering
  
  for (i in 1:length(fixed_vars)) {
  
    # Assert that the default values are actually used  in the simulation.
    # If not, throw an error message so the user can fix it
    # <...> !( fixed_defaults[i] %in% op_chars[fixed_vars[i]] )
    if ( !(  any(fixed_defaults[i] == op_chars[fixed_vars[i]])   )  ) {

      stop(paste('Error: fixed val selected for',fixed_vars[i],'is not in operating characteristics dataframe. select a default that is in the dataframe'  ))
    }
    # Filter to only keep rows with default value of current var
    op_chars <- op_chars %>% filter( eval(as.name(fixed_vars[i])) == fixed_defaults[i] )
    
    # Construct the figure subtitle
    comma<-','
    if (i==length(fixed_vars)) {comma<-''}
    if (i==1) {
      subt <- paste(var_labs[i], fixed_defaults[i],comma)
    }
    else{
      subt <- paste(subt,var_labs[i], fixed_defaults[i],comma)
    }
  }
  
  # Make the plot
  op_chars %>% ggplot(aes(x=eval(as.name(xvar)),
                          y=eval(as.name(yvar)),
                          color=method)) + 
               geom_line() + 
               labs(x=xlab,
                    y=ylab,
                    title=plt_title,
                    subtitle=subt)
            
}

line_plot(op_chars,'u','bias',xlab='Error variance',ylab='Bias',plt_title = 'MSE of ME correction methods, varying error variance',
          fixed_defaults = c(0.3,1,5000,0))

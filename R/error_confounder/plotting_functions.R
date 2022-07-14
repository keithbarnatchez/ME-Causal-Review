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
        plot.caption = ggplot2::element_text(hjust = 0),
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

label_handler <- function(var) {
  # x vars
  if (var=='u') return('Error variance')
  if (var=='a') return('True ATE')
  if (var=='n') return('Sample size')
  if (var=='b') return('Binary outcome')
  if (var=='rho') return('Corr(X,Z)')
  if (var=='psi') return('Corr(X,V)')
  if (var=='ax') return('ax')
  
  # y vars
  if (var=='bias') return('Percent bias')
  if (var=='mse') return('MSE')
  if (var=='ci_cov') return('95% CI coverage rate')
}

title_handler <- function(g) {
  return(1)
}

line_plot <- function(op_chars,xvar,yvar,
                      fixed_vars=c('u','a','n','b','rho','psi','ax'),
                      fixed_defaults=NULL, # c(0.3,1,500,0,0.5,0.05,0.25),
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
  
  # specify defaults to be first row of op_chars if defaults are not provided
  if (is.null(fixed_defaults)) { 
    tempdf <- op_chars %>% select('u','a','n','b','rho','psi','ax') %>%
      unique() 
    fixed_defaults <- tempdf[1,]
  }
  
  # Prune out the var that's being varied
  drop_idx <- (fixed_vars %in% xvar)   
  fixed_vars <- fixed_vars[! drop_idx] ; fixed_defaults <- fixed_defaults[! drop_idx]
  # var_labs filter step was here but changing that
  
  var_labs <- rep(NA,length(fixed_vars)) 
  for (v in 1:length(fixed_vars)) {
    var_labs[v] <- label_handler(fixed_vars[v])
  }
  
  # Fix values at defaults via filtering
  for (i in 1:length(fixed_vars)) { # loop over the variables being kept fixed
  
    # Assert that the default values are actually used  in the simulation.
    # If not, throw an error message so the user can fix it
    # <...> !( fixed_defaults[i] %in% op_chars[fixed_vars[i]] )
    if ( !(  any(as.double(fixed_defaults[i]) == op_chars[fixed_vars[i]])   )  ) {
      stop(paste('Error: fixed val selected for',fixed_vars[i],'is not in operating characteristics dataframe. select a default that is in the dataframe'  ))
    }
    # Filter to only keep rows with default value of current var
    op_chars <- op_chars %>% filter( eval(as.name(fixed_vars[i])) == as.double(fixed_defaults[i]) )
    
    # Construct the figure subtitle
    comma<-', '
    if (i==length(fixed_vars)) {comma<-''}
    if (i==1) {
      subt <- paste(var_labs[i],': ', fixed_defaults[i],comma,sep='')
    }
    else{
      subt <- paste(subt,var_labs[i],': ', fixed_defaults[i],comma,sep='')
    }
  }
  
  # Make the title
  plt_title <- paste(label_handler(yvar),'of ME correction methods')

  # Make the plot
  op_chars %>% ggplot(aes(x=eval(as.name(xvar)),
                          y=eval(as.name(yvar)),
                          color=method)) + 
               geom_line(size=0.75) + geom_point(size=2) +
               labs(x=label_handler(xvar),
                    y=label_handler(yvar),
                    title=plt_title,
                    caption=subt,
                    subtitle=paste('Varying',label_handler(xvar)),
                    color='Method')
  
  # Want to save to current figure 
  # first make the filename (format xvar_yvar_)
  nm_str <- paste(xvar,'_',yvar,sep='')
  for (i in 1:length(fixed_vars)) {
    nm_str <- paste(nm_str,'_',fixed_vars[i],fixed_defaults[i],sep='')
  }
  nm_str <- paste(   gsub('.','',nm_str,fixed=T),'.pdf',sep=''    )
  fullpath <- paste('../../output/figures/',nm_str,sep='')
  ggsave(fullpath,height = 6,width = 9,units = 'in')
}


function (data, columns = 1:ncol(data), color = NULL, alpha = 1, 
          corMethod = "pearson") 
{
  data <- data
  data.choose <- data[columns]
  dn <- data.choose[sapply(data.choose, is.numeric)]
  if (ncol(dn) == 0) {
    stop("All of your variables are factors. Need numeric variables to make scatterplot matrix.")
  }
  if (ncol(dn) < 2) {
    stop("Not enough numeric variables to make a scatter plot matrix")
  }
  a <- uppertriangle(data, columns = columns, color = color, 
                     corMethod = corMethod)
  if (is.null(color)) {
    plot <- scatmat(data, columns = columns, alpha = alpha) + 
      geom_text(data = a, aes_string(label = "r"), colour = "black", family = "serif")
  }
  else {
    plot <- scatmat(data, columns = columns, color = color, 
                    alpha = alpha) + geom_text(data = a, aes_string(label = "r", 
                                                                    color = "colorcolumn")) + labs(color = color)
  }
  factor <- data.choose[sapply(data.choose, is.factor)]
  if (ncol(factor) == 0) {
    return(plot)
  }
  else {
    warning("Factor variables are omitted in plot")
    return(plot)
  }
}

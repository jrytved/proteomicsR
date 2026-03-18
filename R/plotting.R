#' Flip x labels of any plot with 90 degrees
#' @export
flip_xlab <- function(p){
  theme(axis.text.x=element_text(angle=90))
}
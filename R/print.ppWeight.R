#' @export
print.ppWeight<-function(x,...){
  attributes(x) <- NULL
  print(x)
}

#' Visualize estimated membership parameters using ggvis
#' 
#' @param object EstimatedParameters class
#' 
#' @return a figure of estimaged membership parameter via ggvis is generated
#' 
#' @examples 
#' After obtaining EstimatedParameters (typically by \code{getPMSignature}) as Param,
#' visPMSignature2(Param)
#' 
#' @export
setGeneric("visMembership2", function(object) {
  standardGeneric("visMembership2")
})


setMethod("visMembership2",
          signature = c(object = "EstimatedParameters"),
          function(object) {
            
            # function used for add_tooltips
            vis_names <- function(x) {
              if(is.null(x)) return(NULL)
              
              sample <- x[2]
              signature <- x[1]
              intensity <- round(x[4] - x[3], 3)
              paste(sample, signature, intensity, sep="<br/>")
            }
            
            membership_df <- getMembershipValue(object)
            membership_df <- data.frame(sample = rownames(membership_df),  membership_df) %>% gather(signature, intensity, -sample)
            
            
            membership_df %>% ggvis(~sample, ~intensity, fill = ~signature) %>% 
              layer_bars() %>% 
              add_tooltip(vis_names, on = "hover") %>% 
              hide_axis("x")
          
          }
)


            



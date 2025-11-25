import::here(
  "ggplot2", 
  c("element_blank", "element_rect", "element_text", "margin", "%+replace%"),
  .character_only=TRUE
)

# custom ggplot2 theme
theme_custom <- function(){

  # start from a pre-built ggplot2 theme and replace elements we want to change
  ggplot2::theme_linedraw() %+replace%

    ggplot2::theme(
      
      # background elements
      panel.background = element_rect(fill = "white", colour="white"),
      plot.background = element_rect(fill = "white", colour="white"),
      panel.border = element_blank(),

      #title
      plot.title = element_text(
        size = 16,                # set font size
        face = 'bold',            # bold typeface
        hjust = 0,                # left align
        vjust = 1                 # raise slightly
      ),

      #subtitle
      plot.subtitle = element_text(
        size = 15,                # font size
        hjust = 0,
        vjust = 0
      ),
      
      #caption
      plot.caption = element_text(
        size = 13,                # font size
        hjust = 1,                # right align
        vjust = 0
      ), 
      
      #axis title
      axis.title = element_text(
        size = 18,               # font size
      ),
      
      #axis text
      axis.text = element_text(
        size = 18,                # font size
      ),
      
      # legend title
      legend.title = element_text(
        size = 18,                # font size
      ),
      
      # legend text
      legend.text = element_text(
        size = 18,                # font size
      ),
      
      #margin for axis text
      axis.text.x = element_text(
        margin=margin(5, b = 10))
    )
}

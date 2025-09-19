library(ggplot2)

axis_line = element_line(color = "black", linewidth = 1/2.13)
tick_line = element_line(color = "black", linewidth = 1/2.13)
tick_length = unit(3,"pt")
tick_text = element_text(family="sans", face = "bold", color = "black", size = 8)
axis_text = element_text(family="sans", face = "bold", color = "black", size = 10)

geom_noBG = function(){
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
    axis.line = axis_line, axis.ticks=axis_line, axis.ticks.length=tick_length,
    axis.title=axis_text, axis.text=tick_text
  )
}

geom_color = function(cl){
  scale_colour_manual(values = cl)
}

geom_fill = function(cl){
  scale_fill_manual(values = cl)
}

heatmap_color = function(data, midpoint, colorRange,  paletteLength = 20){
  myColor = colorRampPalette(colorRange)(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks = c(seq(min(data), midpoint, length.out=ceiling(paletteLength/2) + 1), 
               seq(midpoint+(max(data)-midpoint)/floor(paletteLength/2), max(data), length.out=floor(paletteLength/2)))
  return(list(colors=myColor, breaks=myBreaks))
}

make_palette=function(x,color.min="#F8E0E0",color.max="darkred"){
  return(colorRampPalette(colors = c(color.min, color.max))(n = length(x)-1))
}

make_palette2=function(x1,x2,x3,color1.min="darkblue",color1.max="lightblue",color.middle="white",color2.min="#F8E0E0",color2.max="darkred"){
  my_palette=c(colorRampPalette(colors = c(color1.min, color1.max))(n = length(x1)-1),
             rep(color.middle,length(x2)+1),
             colorRampPalette(colors = c(color2.min, color2.max))(n = length(x3)-1))
  return(my_palette)
}


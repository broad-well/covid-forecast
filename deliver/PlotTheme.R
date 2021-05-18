library(ggplot2)
library(ggsci)

# Color scheme: ["042a2b","5eb1bf","cdedf6","ef7b45","d84727"]
colors.black <- '#042a2b'
colors.blue <- '#5eb1b4'
colors.cyan <- '#bae6f2'
colors.pink <- '#ef7b45'
colors.red <- '#d84727'

global.theme <- theme_bw(base_family = 'Linux Libertine O', base_size = 13)
  # theme(line = element_line(color = colors.black),
  #       title = element_text(color = colors.black),
  #       panel.border = element_rect(color = '#98daec'),
  #       axis.ticks = element_line(color = '#98daec'),
  #       panel.grid = element_line(color = colors.cyan),
  #       axis.text = element_text(color = '#347783'))
poster.theme <- theme_bw(base_family = 'Linux Libertine O', base_size = 15) +
  theme(plot.background = element_rect('#f1f2f4'))


options(ggplot2.continuous.colour=colors.black)

export <- function (filename, plot = last_plot()) {
  ggsave(filename, plot = plot + 
           scale_color_nejm() +
           scale_fill_nejm(), width = 16/2.4, height = 9/2.4, dpi = 120*2.4)
}

plp <- function(plot) {
  plot + scale_color_nejm() + scale_fill_nejm() + global.theme
}

textless <- function (theme.1) {
  theme.1 + theme(text = element_blank())
}

library(HexSticker)
library(ggplot2)

d <- seq(from = 0, to = 5, by = 0.01)

prior <- median(rlnorm(1E3,1,.5))

truth <- 0.5
truth_shape <- 3.5

df <- tibble(Distance = d,
             Posterior = pweibull(q = d,shape = truth_shape,scale = truth,lower.tail = F),
             Prior = pexp(q = d, rate = prior, lower.tail = F)) %>% 
  gather(Prior,Posterior,key = "Knowledge",value="Estimate")

p <- df %>% ggplot(aes(x = Distance, y = Estimate,linetype = Knowledge)) + geom_line() +
  xlim(0,1) + labs(title = "") + theme_void() + theme(legend.position = "none")
library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Raleway")
## Automatically use showtext to render text for future devices
showtext_auto()

hexSticker::sticker(p, package="rstap", p_size = 7, p_x = 1, p_y = 1.5,
                    s_x=1, s_y=1, s_width=1.3, s_height=1,
                    h_fill = "#f9fcfc", h_color = "#060707",
                    p_color = "black", p_family = "Raleway",
                    filename="Desktop/rstap_hex.png")
    

################################################################################
################### A SCRIPT TO PLAY AROUND WITH IMAGES ########################
################################################################################

pacman::p_load("magick", "hexSticker")

### load global poverty map picture
povimage_png <- image_read("inst/sticker/poverty-map.png")


### create sticker

sticker_obj <- sticker(povimage_png,
                       package = "povmap",
                       p_size = 20,
                       s_x = 1,
                       s_y = 0.85,
                       s_width = 8,
                       s_height = 0.6,
                       h_color = "#89C4F4",
                       h_fill = "#4B77BE",
                       spotlight = TRUE,
                       p_family = "Aller_Lt",
                       filename = "inst/sticker/povmap.png")


usethis::use_logo("inst/sticker/povmap.png")
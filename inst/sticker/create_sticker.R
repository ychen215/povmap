################################################################################
################### A SCRIPT TO PLAY AROUND WITH IMAGES ########################
################################################################################

pacman::p_load("magick", "hexSticker")

### load global poverty map picture
povimage_png <- image_read("inst/sticker/poverty-map.png")
#povimage_png <- image_read("inst/sticker/poverty-map2.jpg")


### create sticker

sticker_obj <-
  sticker(povimage_png,
          package = "povmap",
          p_color = "#FEB845",
          u_color = "white",
          p_size = 80,
          p_y = 1.0,
          s_x = 1,
          s_y = 1,
          s_width = 4,
          s_height = 4,
          white_around_sticker = TRUE,
          h_color = "black",
          dpi = 1000,
          spotlight = TRUE,
          filename = "inst/sticker/sticker_hex.png")

usethis::use_logo("inst/sticker/povmap.png")
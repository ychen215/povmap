
# Generate emdi object with additional indicators; here via function ebp()
emdi_model <- ebp(
  fixed = eqIncome ~ gender + eqsize + cash +
    self_empl + unempl_ben + age_ben + surv_ben + sick_ben +
    dis_ben + rent + fam_allow + house_allow + cap_inv +
    tax_adj, pop_data = eusilcA_pop,
  pop_domains = "district", smp_data = eusilcA_smp,
  smp_domains = "district", threshold = 11064.82,
  transformation = "box.cox", L = 2, MSE = TRUE, B = 2
)


### change the shapefile to an sf object
# shape_austria_dis <- sf::st_as_sf(shape_austria_dis, crs = 4326)
#
# save(shape_austria_dis, file = "inst/shapes/shape_austria_dis.rda")

#### test that the new shapefile works
load_shapeaustria()



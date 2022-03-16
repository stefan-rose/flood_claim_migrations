#Scrpt to run migration models and visualize outputs 
pacman::p_load(tidyverse,tigris,lfe,stargazer,sf,spatialreg,sp,car,
               multiwayvcov,tmap,caret,estimatr,mgcv)

# load in data ----
#natural amenity data supplied by USDA
nat_amenity <- read_csv("~/Desktop/Thesis/Coastal_Zone/natamenf_1_.csv", 
                        col_types = cols(FIPS_Code = col_double()))
#Education data supplied by census bureau
Education <- read_csv("~/Desktop/Thesis/Coastal_Zone/Education.csv", 
                      col_types = cols(FIPS = col_double()))

#US counties supplied by Census TIGER/Lines Shapefiles
county_shp <- counties(cb=T)

#migration claims provided by IRS - with processing by Matt Hauer - Florida State University
#can be downloaded here: https://github.com/mathewhauer/irs-migration-replication
#this data has additional variables added in prior script that includes, FEMA disasters and NFIP flood claims, and demographics by county
mig_clm <- read_csv("~/Desktop/Thesis/IRS_mig_data/mig_claims_4_22.csv")

mig_clm_shp <- left_join(mig_clm,county_shp,by=c("fips"="GEOID")) %>% st_as_sf()
mig_clm_shp <- mig_clm_shp %>%  select(-STATEFP.y:-AWATER.y)

mig_clm_shp <- mig_clm_shp %>% left_join(nat_amenity,by=c("fips"="FIPS_Code"))

#redefine Great Lakes States as non-coastal
mig_clm_shp <- mig_clm_shp %>% mutate(coastal = if_else(state_name == "Michigan"|
                                                          state_name== "Minnesota"|
                                                          state_name=="Wisconsin"|
                                                          state_name=="Illinois"|
                                                          state_name=="Ohio"|
                                                          state_name=="Indiana"|
                                                          state_name=="Pennsylvania",0,coastal)) %>% 
  mutate(coastal = if_else(fips == 36029|fips==36045|fips==36055|fips==36063|fips==36073|
                             fips==36075|fips==36117,0,coastal))

mig_fil <- mig_clm_shp %>% ungroup() %>% group_by(fips) %>% 
  summarise(flood_disaster = sum(flood_disaster),
            other_disaster = sum(other_disaster),
            out_flow = mean(out_flow, na.rm=T),
            amenity_rank = last(nat_amenity_rank),
            population = last(population))

#vizualize maps of disasters per capita 
tm_shape(county_shp, projection = proj) + tm_polygons(col="lightgrey",alpha = 0.9,border.col = "white",border.alpha = 0.5) +
  tm_shape(mig_fil,is.master = TRUE, projection = proj) + 
  tm_polygons("disasters_pc",style="quantile", border.col=NULL,
              alpha=0.8,palette="PuRd",title="FEMA Disasters Declared") +
  tm_layout(frame = FALSE)

png("disasters_pc.png",pointsize = 11,height = 900,width = 1150,res = 150)
dev.off()


# base model specifications -----
M1 <- lm(log10(out_flow)~flood_disaster+other_disaster+total_claims+coastal+HPI+Unemployment_Rate+Jan_tmin+as.factor(year)+
           as.factor(state_name), data=mig_clm_shp)

M2 <- lm(log(out_flow)~flood_disaster+other_disaster+total_claims+coastal+nat_amenity_rank+
           HPI+Unemployment_Rate+Jan_tmin+log(pop_density),
         data=mig_clm_shp, na.action = NULL)

summary(M4)


#fixed effects by state and year
mig.fe_lfe <- felm(log(out_flow)~flood_disaster+other_disaster+total_claims+coastal+nat_amenity_rank+
                     HPI+Unemployment_Rate+Jan_tmin+log(pop_density)|
                     state_name+year|0|state_name,
                   data = mig_clm_shp,na.action = NULL)

#fe year
mig.fe_year <- felm(log(out_flow)~flood_disaster+other_disaster+total_claims+coastal+nat_amenity_rank+
                      HPI+Unemployment_Rate+Jan_tmin+log(pop_density)|
                      year|0|0,
                    data = mig_clm_shp,na.action = NULL)
#fe state
mig.fe_state <- felm(log(out_flow)~flood_disaster+other_disaster+total_claims+coastal+nat_amenity_rank+
                       HPI+Unemployment_Rate+Jan_tmin+log(pop_density)+year|
                       state_name|0|0,
                     data = mig_clm_shp,na.action = NULL)
#fe coastal zone lag added
mig.fe_coastal <- felm(log10(out_flow)~flood_disaster+other_disaster+total_claims+HPI+Unemployment_Rate+Jan_tmin+lag_s_out_flow|
                         coastal|0|0,
                       data = mig_clm_shp,na.action = NULL)
#fe year manual lag added
mig.fe_lag <- felm(log10(out_flow)~flood_disaster+other_disaster+total_claims+HPI+Unemployment_Rate+Jan_tmin+lag_s_out_flow|
                     year|0|0,
                   data = mig_clm_shp,na.action = NULL)
mig.ols_trend <- felm(out_flow~flood_disaster+other_disaster+total_claims+HPI+Unemployment_Rate+Jan_tmin+year+(year^4)|
                        0|0|0,
                      data = mig_clm_shp,na.action = NULL)
#no fe, standard ols
mig.ols_lfe <- felm(log(out_flow)~flood_disaster+other_disaster+total_claims+coastal+nat_amenity_rank+
                      HPI+Unemployment_Rate+log(pop_density)+year
                    |0|0|0,
                    data = mig_clm_shp,na.action = NULL)
summary(mig.fe_lfe)
#create output tables with stargazer
multiply.100 <- function(x) (x*100) #function to make coefficient and error interpretations easier in outputs

stargazer(M1,mig.ols_lfe,mig.fe_year,mig.fe_state,mig.fe_coastal,mig.fe_lfe,type = "html",
          model.names = FALSE,
          add.lines = list(c("Year FE","-","-","Yes","-","Yes","Yes"),
                           c("Region FE","-","-","-","Yes","Yes","Yes"),
                           c("Corrected Spatial Autocorrelation","-","-","-","-","-","Yes")),
          out = "reg_table4.htm")


stargazer(mig.ols_lfe,mig.fe_year,mig.fe_state,mig.fe_lfe,type = "html",
          title = "Modeling Migration Outflows with Natural Disaster Indicators",
          model.names = FALSE,omit.stat = c("ser","f"),
          apply.coef = multiply.100, apply.se = multiply.100,
          add.lines = list(c("Year FE","-","Yes","-","Yes"),
                           c("Region FE","-","-","Yes","Yes")),
          out = "clusters_table2.htm")

stargazer(mig.ols_lfe,mig.fe_year,mig.fe_state,mig.fe_lfe,type = "html",
          title = "Modeling Migration Outflows with Natural Disaster Indicators",
          model.names = FALSE,omit.stat = c("ser","f"),
          apply.coef = multiply.100, apply.se = multiply.100,
          add.lines = list(c("Year FE","-","Yes","-","Yes"),
                           c("Cluster FE","-","-","Yes","Yes")),
          out = "clus_fe_table1.htm")


#dot whisker plots to visualize coefficients between models
dwplot(list(M1,mig.ols_lfe,mig.fe_year,mig.fe_state,mig.fe_coastal,mig.fe_lfe)) + 
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Modeling Disaster Mobility") +
  theme(plot.title = element_text(face="bold"),
        legend.position = c(0.80, 0.65),
        legend.justification = c(0, 0), 
        legend.background = element_rect(colour="grey80"),
        legend.title = element_blank())

dwplot(list(mig.ols_lfe,mig.fe_year,mig.fe_state,mig.fe_lfe)) + 
  theme_bw() + xlab("Coefficient Estimate") + ylab("") +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("Modeling Disaster Mobility (Clustered FE)") +
  theme(plot.title = element_text(face="bold"),
        legend.position = c(0.80, 0.65),
        legend.justification = c(0, 0), 
        legend.background = element_rect(colour="grey80"),
        legend.title = element_blank())

########Spatial Models#####
#manual demean year for year fe in spatial model
mig_year_demeaned <- with(mig_clm_shp,
                          data.frame(out_flow = log(out_flow) - ave(log(out_flow), year),#-ave(log(out_flow),state_name),
                                     net_flow = net_flow - ave(net_flow,year),#-ave(net_flow,state_name),
                                     in_flow = log(in_flow) - ave(log(in_flow),year),#-ave(log(in_flow),state_name),
                                     total_movement = log(total_movement) - ave(log(total_movement),year),#-ave(log(total_movement),state_name),
                                     # disaster = disaster - ave(disaster,year)-ave(disaster,nat_amenity_rank),
                                     flood_disaster = flood_disaster - ave(flood_disaster, year),#-ave(flood_disaster,state_name),#-ave(flood_disaster,state_name),
                                     other_disaster = other_disaster - ave(other_disaster, year),#-ave(other_disaster,state_name),#-ave(other_disaster,state_name),
                                     total_claims = total_claims - ave(total_claims, year),#-ave(total_claims,state_name),#-ave(total_claims,state_name),
                                     coastal = coastal - ave(coastal, year),#-ave(coastal,state_name),
                                     amenity_rank = nat_amenity_rank-ave(nat_amenity_rank,year),#-ave(nat_amenity_rank,state_name),
                                     HPI = HPI - ave(HPI, year),#-ave(HPI,state_name),#-ave(HPI,state_name),
                                     Unemployment_Rate = Unemployment_Rate - ave(Unemployment_Rate, year),#-ave(Unemployment_Rate,state_name),#-ave(Unemployment_Rate,state_name),
                                     Jan_tmin = Jan_tmin - ave(Jan_tmin, year),#-ave(Jan_tmin,state_name),
                                     pop_density = log(pop_density) - ave(log(pop_density), year),#-ave(log(pop_density),state_name),#-ave(Jan_tmin,state_name),
                                     year=year,
                                     geometry=geometry))

# add spatial data and create spatial neighborhood structure ---------
mig_clm_shp <- st_as_sf(mig_clm_shp)

mig.nb <- spdep::poly2nb(mig_clm_sp,queen = TRUE)
mig.listw <- spdep::nb2listw(mig.nb,style="W",zero.policy=TRUE)

#spatial two stage least squares with manual state and year fe 
stsls.mig <- stsls(out_flow~flood_disaster+other_disaster+total_claims+
                     coastal+HPI+Unemployment_Rate+Jan_tmin, data = mig_year_demeaned,listw = mig.listw,
                   na.action = NULL,robust = TRUE)
summary(stsls.mig)

#generalized spatial two stage least squares with manual year fe
gstsls.mig <- spatialreg::gstsls(out_flow~flood_disaster+other_disaster+total_claims+coastal+amenity_rank+
                                   HPI+Unemployment_Rate+Jan_tmin+pop_density
                                 , data = mig_year_demeaned,listw = mig.listw,
                                 na.action = NULL,robust = TRUE)
summary(gstsls.mig)

#test for autocorrelation 
spdep::moran.mc(gstsls.mig$resid,mig.listw,nsim = 99,zero.policy = T)




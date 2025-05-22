# Index for body weight maintenance was calculated as described in
# "Variation in extracellular matrix genes is associated with 
#  weight regain after weight loss in a sex-specific manner"

#-------------------------------------------------------------------------------
# Define functions to calculate indexes to base responders and non responders

index_cid4_cid2_fun <- function(x, y){(x[y == "CID_4"] - x[y == "CID_2"])}

index_wms_fun <- function(x, y){(x[y == "CID_4"] - x[y == "CID_2"])/
    (x[y == "CID_1"] - x[y == "CID_2"])}


#-------------------------------------------------------------------------------
# Define functions to assign responders and non_responders based on index
#-------------------------------------------------------------------------------
resp_5tile_fun <- function(x){case_match(ntile(x, 5), 
                                         1 ~ "Resp", 2 ~ "Resp",
                                         3 ~ NA,
                                         4 ~ "NonResp", 5 ~ "NonResp")}

resp_3tile_fun <- function(x){case_match(ntile(x, 5), 
                                         1 ~ "Resp", 
                                         2 ~ NA,
                                         3 ~ "NonResp")}

resp_morethan0_fun <- function(x){ifelse(x <= 0, "Resp", "NonResp")}


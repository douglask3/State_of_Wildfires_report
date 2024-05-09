graphics.off()
library(terra)
library(modi)

regions = c('Greece', 'Canada', 'NW_Amazon')
date_tests = c('2023-08', '2023-06', '2023-09')
out_dir = "outputs/obs_time_series/"

for_region <- function(region, date_test) {
    month_number = as.numeric(substr(date_test, 6, 7))
    burnt_area_data = paste0("data/data/driving_data/", region, "/isimp3a/obsclim/GSWP3-W5E5/period_2000_2019/burnt_area-2000-2023.nc")
    date_test = '2023-06'
    burnt_area = rast(burnt_area_data)
    date_test = substr(time(burnt_area), 1, 7) == date_test
    #burnt_area_event = burnt_area[[date_test]]
    gridArea = cellSize(burnt_area[[1]])
    vArea = values(gridArea)
    mean_95 <- function(i) {
        #weighted.quantile(burnt_area[[i]][], gridArea[[i]], prob = 0.95, plot = FALSE)
        print(i)
        vr = values(burnt_area[[i]])
        val = weighted.quantile(vr, vArea, prob = 0.95, plot = FALSE)
    
        test = vr >= val & !is.na(vr)
        out = sum(vr[test] * vArea[test])/sum(vArea[test])
        return(out)
    }

    burnt_area_tot = sapply(1:nlyr(burnt_area), mean_95)
    names(burnt_area_tot) = time(burnt_area)
    out_dir_region = paste0(out_dir, '/', region, '/')    
    try(dir.create(out_dir_region, recursive = TRUE))
    out_file = paste0(out_dir_region, 'burnt_area_95_monthly.csv')
    write.csv(burnt_area_tot, file = out_file)
    
    burnt_area_event = burnt_area_tot[date_test]
    burnt_area_tot = burnt_area_tot[sort(unlist(lapply((month_number-1):(month_number+1), function(i) seq(i, nlyr(burnt_area), by = 12))))]
    percentile = mean(burnt_area_tot <= burnt_area_event)
    
    out_file = paste0(out_dir_region, 'outs.Rd')
    save(burnt_area_tot, percentile, burnt_area_event, file = out_file)
}

mapply(for_region, regions, date_tests)

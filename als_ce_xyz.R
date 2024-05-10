require(getopt, quietly = TRUE)

spec = matrix(c(
  'help'       , 'h', 0, "logical",
  'input'      , 'i', 1, "character", ## -- path with ALS files
  'output'     , 'o', 1, "character", ## -- path to write structural complexity metrics
  'plots_path' , 'f', 1, "character", ## -- input plot locations (geospatial points file, e.g. GEDI footprint coordinates)
  'n_plots'    , 'n', 1, "integer", ## -- if no input plot locations are provided, how many random samples to draw?
  'plot_size'  , 's', 1, "integer", ## -- plot diameter, in point cloud units (e.g. 25m for GEDI footprints)
  'gridded'    , 'g', 0, "logical", ## -- sample on a regular grid instead of randomly
  'las_epsg'   , 'e', 1, "integer"  ## -- EPSG code of ALS point clouds (if not in LASheader)
), byrow=TRUE, ncol=4)

opt = getopt(spec)

if ( !is.null(opt$help) || is.null(opt$input) || is.null(opt$output) ) {
  cat(getopt(spec, usage=TRUE))
  q(save='no',status=1)
}

PLOT_SIZE = 25
if(!is.null(opt$plot_size)) PLOT_SIZE = opt$plot_size

N_PLOTS = 100
if(!is.null(opt$n_plots)) N_PLOTS = opt$n_plots

GRIDDED = FALSE
if(!is.null(opt$gridded)) GRIDDED = TRUE

PLOTS_PATH = NULL
if(!is.null(opt$plots_path)) PLOTS_PATH = opt$plots_path

LAS_EPSG = NA
if(!is.null(opt$las_epsg)) LAS_EPSG = opt$las_epsg

IN_PATH = opt$input
OUT_PATH = opt$output

## -- load libraries
require(future, quietly = TRUE)
require(future.apply, quietly =  TRUE)
require(lidR, quietly = TRUE)
require(magrittr, quietly = TRUE)
require(TreeLS, quietly = TRUE) # remotes::install_github('tiagodc/TreeLS')
require(data.table, quietly = TRUE)
require(wk, quietly = TRUE)
require(sf, quietly = TRUE)
require(terra, quietly = TRUE)

## -- 3D Canopy Entropy (Liu et al. 2022)
avg_point_dist = function(las, h=1){
  las@data$MinDist = nabor::knn(TreeLS:::las2xyz(las), k=2)$nn.dist[,2]
  z_pts = las$Z - min(las$Z)
  las@data$layer = as.integer(z_pts / h)
  avg_dst = las@data[order(layer),.(dist = mean(MinDist)),by=layer]
  return(avg_dst)
}

adaptative_resample = function(las, h=1, alpha = 0.05){
  avg_dst = avg_point_dist(las, h)
  
  noise_tol = mean(avg_dst$dist) + 2*sd(avg_dst$dist)
  avg_dst = avg_dst[dist < noise_tol,]
  max_dst = max(avg_dst$dist)
  
  for (ratio in seq(2,5,.1)){
    vox = ratio * max_dst
    vlas = tlsSample(las, smp.voxelize(vox))
    avg_dst = avg_point_dist(vlas, h)
    mkt = trend::mk.test(avg_dst$dist)
    if(mkt$p.value > alpha) break
  }
  
  return(list(las=vlas,ratio=ratio,vox=vox,dist=avg_dst,mk.test=mkt))
}

canopy_entropy = function(las, bw=.2, grid_size=.1){
  bounds = apply(las@data, 2, range)
  ce = c()
  for ( ij in list(c('X','Y'), c('X','Z'), c('Y','Z')) ){
    i = ij[1]
    j = ij[2]
    plane = las@data[,..ij]
    
    ni = 2 + (diff(bounds[,i]) + 8*bw) %/% grid_size
    nj = 2 + (diff(bounds[,j]) + 8*bw) %/% grid_size
    lims = c(bounds[,ij]) + (bw*4 +(grid_size/2))*c(-1,1,-1,1)
    
    den = ks::kde(plane, h=bw, gridsize = c(ni,nj), xmin = lims[c(1,3)], xmax = lims[c(2,4)])
    den = den$estimate[den$estimate > 0]
    entropy = -1 * sum(den*log(den)*grid_size*grid_size)
    
    ce[paste(ij,collapse = '')] = entropy
  }

  ce['XYZ'] = sqrt(sum(ce^2))
  df = ce %>% t %>% as.data.frame
  return(df)
}

get_entropy = function(las, h=1){
    las = filter_poi(las, Height >= h)
    if(is.empty(las)) return(NULL)
    las@data$Z = las@data$Height
    res = adaptative_resample(las)
    
    ent = canopy_entropy(res$las)
    ent$p = res$mk.test$p.value
    ent$vox = res$vox

    return(ent)
}

## -- pipeline
las_height = function(las, dtm_res=0.5){
    dtm = rasterize_terrain(las, res=dtm_res, algorithm = knnidw())
    las = normalize_height(las, dtm)
    las = add_lasattribute(las, las@data$Z, 'Height', "Height above ground")
    las = unnormalize_height(las)    
    return(las)
}

get_complexity = function(las){
    if( !('Height' %in% names(las@data)) ){
        las = las_height(las)
    } 
    ce = get_entropy(las)
    return(ce)
}

clip_and_process = function(pt, ctg, l=PLOT_SIZE){
    x = as.double(pt %>% st_coordinates)[1]
    y = as.double(pt %>% st_coordinates)[2]
    
    buff = pt %>% st_buffer(l/2)
    inter = st_intersects(buff, ctg@data)[[1]]
    files = ctg@data[inter,]$filename
   
    filt = glue::glue("-keep_circle {x} {y} {l/2}")
    filt = paste(filt, opt_filter(ctg))
    cols = opt_select(ctg)
    
    if(length(files) == 0) return(NULL)    
    
    las = readLAS(files, select = cols, filter = filt)
    
    if(is.empty(las)) return(NULL)
    if(nrow(las@data) < 100) return(NULL)
    if(!any(las$Classification == 2)) return(NULL)
    
    comp = tryCatch(get_complexity(las), error=function(e) NULL)
    return(comp)
    
}

## -- run
n_cores = parallel::detectCores()
n_workers = as.integer(n_cores / 4)
set_lidr_threads(1)
plan(multicore, workers = n_workers)

cat(glue::glue('\n## -- calculating complexity metrics for {IN_PATH}\n'))

ctg = readLAScatalog(IN_PATH)
opt_chunk_buffer(ctg) = 0 
opt_stop_early(ctg) = FALSE

if(is.na(st_crs(ctg))){
    st_crs(ctg) = LAS_EPSG
}                    
                    
lhd = readLASheader(ctg$filename[1])
h_byte = which(names(lhd@VLR$Extra_Bytes$`Extra Bytes Description`) == 'Height')
opt_select(ctg) = paste0('xyzc', h_byte)

buff = ctg@data %>% st_union %>% st_buffer(-PLOT_SIZE/2)

geo_index = NULL
if(is.null(PLOTS_PATH)){
    if(GRIDDED){
        sample_pts = st_sample(buff, size=as.double(st_area(buff)) %/% (PLOT_SIZE^2), type = 'regular')
        if(N_PLOTS > 0) sample_pts = sample_pts %>% sample(N_PLOTS)
    }else{
        sample_pts = st_sample(buff, size=N_PLOTS, type ='random')
    }
    st_crs(sample_pts) = st_crs(ctg)
}else{
    sample_pts = st_read(PLOTS_PATH)
    geo_index = sample_pts$index    
    sample_pts = st_transform(sample_pts, st_crs(ctg))$geom
}

cat(glue::glue('\n## -- opening {n_workers} parallel processes to process {length(sample_pts)} plots\n'))
complexity = future_lapply(sample_pts, clip_and_process, ctg=ctg, l=PLOT_SIZE, future.seed=TRUE)

keep = !sapply(complexity, function(x) is.null(x) || nrow(x) == 0)
cpx = do.call(rbind, complexity)
cpx$geometry = sample_pts[keep]

cpx = st_as_sf(cpx)

if(!is.null(geo_index)){
    cpx$index = geo_index[keep]
}
                   
st_write(cpx, OUT_PATH)
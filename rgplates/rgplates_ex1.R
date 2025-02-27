
# Install 'rgplates' in R:
install.packages("geojsonsf")
install.packages('rgplates', repos = c('https://gplates.r-universe.dev', 'https://cloud.r-project.org'))

library(geojsonsf)
library(rgplates)

# Help
?rgplates

# ?reconstrucdt
# "PALEOMAP" (Scotese, 2016) for coastlines only (0-1100 Ma).

# Takes ~10 seconds each
paleomap_000Ma = reconstruct(x="coastlines", age=00, model="PALEOMAP")
paleomap_080Ma = reconstruct(x="coastlines", age=80, model="PALEOMAP")
paleomap_160Ma = reconstruct(x="coastlines", age=160, model="PALEOMAP")

# "sf"         "data.frame"
class(tmp_map)

# "{sf} spatial objects are simply data.frames with a sticky geometry column"
# https://ucd-cws.github.io/CABW2020_R_training/m2_3_using_sf.html
summary(paleomap_160Ma)

# Plots with lots of lines work better inside of PDFs

pdffn = "paleomaps_v1.pdf"
pdf(file=pdffn, width=6, height=6)

plot(paleomap_000Ma); title("Scotese Paleomap at 0 Ma")

plot(paleomap_080Ma); title("Scotese Paleomap at 80 Ma")

plot(paleomap_160Ma); title("Scotese Paleomap at 160 Ma")

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


# Check projection (CRS=coordinate reference system)
# https://ucd-cws.github.io/CABW2020_R_training/m2_3_using_sf.html
# https://r.geocompx.org/spatial-class.html#crs-intro
st_crs(paleomap_160Ma)

# Let's change to Azimuthal Equal Area projection
# view it:
# get the code, google "EPSG Lambert Azimuthal Equal Area"
# answer: EPSG:9820

st_crs("ESRI:54030")
st_crs("EPSG:9820")



#######################################################
# Projections are crazy complex, it takes awhile to figure out which projections
# (A) are useful/pretty
# (B) work with your particular software package
#
# ...some experiments below. Not worth a huge amount of time, the real project is to 
# extract & consistently label the polygons.
#######################################################

#######################################################
# Kinda worked
#######################################################

# Reproject to a different CRS (e.g., EPSG 3395 - Albers Equal Area)
paleomap_000Ma_AlbertsEqualArea = st_transform(x=paleomap_000Ma, crs="EPSG:3395")
paleomap_080Ma_AlbertsEqualArea = st_transform(x=paleomap_080Ma, crs="EPSG:3395")
paleomap_160Ma_AlbertsEqualArea = st_transform(x=paleomap_160Ma, crs="EPSG:3395")


pdffn = "paleomaps_v1_AlbertsEqualArea.pdf"
pdf(file=pdffn, width=6, height=6)

plot(paleomap_000Ma_AlbertsEqualArea); title("Scotese Paleomap at 0 Ma, AlbertsEqualArea projection")

plot(paleomap_080Ma_AlbertsEqualArea); title("Scotese Paleomap at 80 Ma, AlbertsEqualArea projection")

plot(paleomap_160Ma_AlbertsEqualArea); title("Scotese Paleomap at 160 Ma, AlbertsEqualArea projection")

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)




# WGS 84 / Antarctic Polar Stereographic - EPSG:3031
paleomap_000Ma_Antarctic_Polar = st_transform(x=paleomap_000Ma, crs="EPSG:3031")
paleomap_080Ma_Antarctic_Polar = st_transform(x=paleomap_080Ma, crs="EPSG:3031")
paleomap_160Ma_Antarctic_Polar = st_transform(x=paleomap_160Ma, crs="EPSG:3031")


pdffn = "paleomaps_v1_Antarctic_Polar.pdf"
pdf(file=pdffn, width=6, height=6)

plot(paleomap_000Ma_Antarctic_Polar); title("Scotese Paleomap at 0 Ma, Antarctic_Polar projection")

plot(paleomap_080Ma_Antarctic_Polar); title("Scotese Paleomap at 80 Ma, Antarctic_Polar projection")

plot(paleomap_160Ma_Antarctic_Polar); title("Scotese Paleomap at 160 Ma, Antarctic_Polar projection")

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)








#######################################################
# DIDN'T WORK
#######################################################

paleomap_000Ma_projected = st_transform(x=paleomap_000Ma, pipeline ="+proj=natearth")
paleomap_080Ma_projected = st_transform(x=paleomap_080Ma, pipeline ="+proj=natearth")
paleomap_160Ma_projected = st_transform(x=paleomap_160Ma, pipeline ="+proj=natearth")

pdffn = "paleomaps_v1_projected.pdf"
pdf(file=pdffn, width=6, height=6)

plot(paleomap_000Ma_projected); title("Scotese Paleomap at 0 Ma, projected")

plot(paleomap_000Ma_projected); title("Scotese Paleomap at 80 Ma, projected")

plot(paleomap_000Ma_projected); title("Scotese Paleomap at 160 Ma, projected")

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)




# sf projections available are those in PROJ:
# https://proj.org/en/9.3/operations/projections/index.html
# Azimuthal Equidistant - ESRI:54032
paleomap_000Ma_projected = st_transform(x=paleomap_000Ma, crs="EPSG:27200")
paleomap_080Ma_projected = st_transform(x=paleomap_080Ma, crs="EPSG:27200")
paleomap_160Ma_projected = st_transform(x=paleomap_160Ma, crs="EPSG:27200")

pdffn = "paleomaps_v1_AlbertsEqualArea.pdf"
pdf(file=pdffn, width=6, height=6)

plot(paleomap_000Ma_projected); title("Scotese Paleomap at 0 Ma, projected")

plot(paleomap_000Ma_projected); title("Scotese Paleomap at 80 Ma, projected")

plot(paleomap_000Ma_projected); title("Scotese Paleomap at 160 Ma, projected")

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)


# Mollweide - ESRI:54009
# Reproject to a different CRS (e.g., EPSG 3395 - Albers Equal Area)
paleomap_000Ma_Mollweide = st_transform(x=paleomap_000Ma, crs="EPSG:3395")
paleomap_080Ma_Mollweide = st_transform(x=paleomap_080Ma, crs="EPSG:3395")
paleomap_160Ma_Mollweide = st_transform(x=paleomap_160Ma, crs="EPSG:3395")


pdffn = "paleomaps_v1_Mollweide.pdf"
pdf(file=pdffn, width=6, height=6)

plot(paleomap_000Ma_Mollweide); title("Scotese Paleomap at 0 Ma, AzEqArea projection")

plot(paleomap_080Ma_Mollweide); title("Scotese Paleomap at 80 Ma, AzEqArea projection")

plot(paleomap_160Ma_Mollweide); title("Scotese Paleomap at 160 Ma, AzEqArea projection")

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)



# Gall-Peters EPSG:54016
paleomap_000Ma_Gall = st_transform(x=paleomap_000Ma, crs="EPSG:54016")
paleomap_080Ma_Gall = st_transform(x=paleomap_080Ma, crs="EPSG:54016")
paleomap_160Ma_Gall = st_transform(x=paleomap_160Ma, crs="EPSG:54016")


pdffn = "paleomaps_v1_Gall.pdf"
pdf(file=pdffn, width=6, height=6)

plot(paleomap_000Ma_Gall); title("Scotese Paleomap at 0 Ma, Gall projection")

plot(paleomap_080Ma_Gall); title("Scotese Paleomap at 80 Ma, Gall projection")

plot(paleomap_160Ma_Gall); title("Scotese Paleomap at 160 Ma, Gall projection")

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)




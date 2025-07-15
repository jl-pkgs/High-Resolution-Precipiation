# using GeoStats
using Meshes
# using Meshed: LatLon
using GeoStatsFunctions, Distances
using GeoStatsFunctions: fit
using GeoStatsModels
using SpatRasters
import SpatRasters: st_dims
using RTableTools
using GeoTables
# using GLMakie, MakieLayers

SK = GeoStatsModels.SimpleKriging
OK = GeoStatsModels.OrdinaryKriging
UK = GeoStatsModels.UniversalKriging
DK = GeoStatsModels.UniversalKriging

sf_point(x, y) = Meshes.Point(Meshes.LatLon(y, x)) # note `lat` at first

function st_dims(domain::RegularGrid)
  p0 = domain.origin.coords
  y0, x0 = p0.lat.val, p0.lon.val

  ny, nx = size(domain.topology)
  celly, cellx = domain.spacing
  celly = celly.val
  cellx = cellx.val

  lon = (x0 + cellx / 2):cellx:(x0 + cellx * nx) #|> collect
  lat = (y0 + celly / 2):celly:(y0 + celly * ny) #|> collect
  lon, lat#[end:-1:1]
end


b = bbox(108.5, 31.5, 112 - 0.5, 33.5)
cellsize = 0.02
domain = RegularGrid(sf_point(b.xmin, b.ymin), sf_point(b.xmax, b.ymax), (1, 1) .* cellsize)
p = centroid(domain, 1)


d = fread("data/prcp_st174_shiyan.csv")
data = georef(d[:, [1, 2, 4]], (:lon, :lat))


var = :prcp
γ_emp = EmpiricalVariogram(data, var, maxlag=200_000.0)
γ = GeoStatsFunctions.fit(GaussianVariogram, γ_emp) # `maxlag` in the unit of m
ok = GeoStatsModels.fit(OK(γ), data)
uk = GeoStatsModels.fit(UK(γ, 2, 2), data)
# exponents(2, 2)
point = true
prob = false

# @time r = GeoStatsModels.fitpredict(ok.model, data, domain; point, prob, neighbors=false, distance=Haversine())
r = GeoStatsModels.fitpredictfull(uk.model, data, domain, point, prob)


ra = georef(r, domain)
A = reshape(ra.prcp, size(domain))' |> collect
lon, lat = st_dims(domain)


begin
  colors = amwg256
  colorrange = (0, 40)

  fig = Figure(; size=(800, 400))
  # ax = Axis(fig[1, 1])
  ax, plt = imagesc!(fig[1, 1], lon, lat, A; colors, colorrange, force_show_legend=true)
  scatter!(ax, d.lon, d.lat; color=d.prcp, strokecolor=:white, strokewidth=1, colormap=colors, colorrange)
  fig

  save("ex-UK.png", fig)
end

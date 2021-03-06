## Land-based Station: ##
* **Quality Controlled Local Climatological Data (QCLCD)**

Experiencing a quality controlled, the QCLCD data technically is from Jan 2005 to present. But
since data from 2005 to 2007 is still under processing of quality controlled, the data is from
QCLCD200705 to QCLCD201410, and is consist of hourly, daily, and monthly summaries for
approximately 1,600 U.S. locations. Totally 90 months. Each hourly data (for each month) is about
300Mb, and daily data is about 4Mb, and monthly data is about 200Kb. For hourly data, there are
24 columns, such as DryBulb temperature, WetBulb temperature, DewPoints temperature, RelativeHumidity,
Windspeed, so forth.

* **Unedited Local Climatological Data (ULCD)**

Unedited local data is from July 1996 to July 2007. For over 700 U.S. locations. Each hourly data
(for each month) is about 45 Mb. If we would like to study this LCD data, hourly data is more
reasonable with respect to the size, it is large. And unedited LCD can be merged to QCLCD, which
means we can have data from Junly 1996 to Oct 2014 which is 220 months. But the number of stations
is very small.

* Integrated Surface Database (ISD)

* Global Historical Climatology Network (GHCN)


## Satellite Data: ##

* **Geostationary IR Channel Brightness Temperature - GridSat B1**

There are three variables in this satellite data, infrared window observation for brightness temperature,
infrared water vapor observation, and visible channel observation. Spatial resolution is about 8km
(which is about 0.07 degree in lon and lat) from 70S to 70N. For 1980 to present, every 3 hour observations.

## Reanalysis Data: ##

Reanalysis datasets are created by assimilating ("inputting") climate observations using the same
climate model throughout the entire reanalysis period in order to reduce the affects of modeling
changes on climate statistics. Observations are from many different sources including ships,
satellites, ground stations, RAOBS, and radar.

Reanalysis of past weather data presents a clear picture of past weather, independent of the many
varieties of instruments used to take measurements over the years. Through a variety of methods,
observations from various instruments are added together onto a regularly spaced grid of data.
Placing all instrument observations onto a regularly spaced grid makes comparing the actual
observations with other gridded datasets easier. In addition to putting observations onto a grid,
reanalysis also holds the gridding model constant—it doesn't change the programming—keeping the
historical record uninfluenced by artificial factors. Reanalysis helps ensure a level playing field
for all instruments throughout the historical record.

### Key Strengths: ###
* Global data sets, consistent spatial and temporal resolution over 3 or more decades, hundreds of
variables available; model resolution and biases have steadily improved
* Reanalyses incorporate millions of observations into a stable data assimilation system that would
be nearly impossible for an individual to collect and analyze separately, enabling a number of
climate processes to be studied
* Reanalysis data sets are relatively straightforward to handle from a processing standpoint
(although file sizes can be very large)

### Key Limitations: ###

* Reanalysis data sets should not be equated with "observations" or "reality"
* The changing mix of observations, and biases in observations and models, can introduce spurious
variability and trends into reanalysis output
* Observational constraints, and therefore reanalysis reliability, can considerably vary depending
on the location, time period, and variable considered

### Dataset ###

* **NCEP North American Regional Reanalysis: NARR**

NCEP's high resolution combined model and assimilated dataset. It covers 1979 to near present and is
provided 8-times daily, daily and monthly on a Northern Hemisphere Lambert Conformal Conic grid for
all variables. Current total of 29.4 Tbytes for 8xDaily. The grid resolution is 349x277 which is
approximately 0.3 degrees (32km) resolution at the lowest latitude. Three different type of responses
in the dataset.

1. Air Temperature, Humidity, Zonal and Meridional Wind at 29 pressure levels.
2. Air Temperature at 2 meter, at surface, Accumulated total evaporation at surface, Moisture
Availability at surface, Precipitation rate at surface, Pressure at surface and so forth.
3. Soil Temperature at 5 depth levels.
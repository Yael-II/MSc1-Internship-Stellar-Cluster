"""
* SCOPE : Sky Coordinates for Observations Python Estimator
* Version 1 - 2024-02-09
@ Yaël Moussouni
@ Unistra, P&E, MSc1-MoFP
@ Observatory of Strasbourg (Intern)
"""
import astropy.time as time
import astropy.coordinates as coord
import astropy.units as u

# * Params
eps = coord.Angle("0d0m1s")
# * Location
obs_name = "Strasbourg Observatory"
obs_lat = coord.Latitude("""48d34m59s""")  # @ 48°34'59" N
obs_lon = coord.Longitude("""07d46m07s""") # @ 07°46'07" E

# * User input
print("")
print("\033[94m"+"Info: Location is manually set."+"\033[0m")
print("")
min_alt = coord.Angle(input("\033[32m"+"Minimum altitude above horizon (deg): "+"\033[0m")+"d")
min_dec = coord.Angle(input("\033[32m"+"Minimum declination around the north (deg): "+"\033[0m")+"d")
window_east = coord.Angle(input("\033[32m"+"Observation window before crossing meridian (h): "+"\033[0m")+" hours")
window_west = coord.Angle(input("\033[32m"+"Observation window after crossing meridian (h): "+"\033[0m")+" hours")
obs_date = input("\033[32m"+"Date of observation (YYYY-MM-DD format): "+"\033[0m")
obs_utc = input("\033[32m"+"Time of observation (hh:mm or hh:mm:ss formats - UTC): "+"\033[0m")
obs_duration = coord.Angle(input("\033[32m"+"Observation duration (h): "+"\033[0m")+" hours")
obs_sky_rotation = obs_duration * (1*u.sday/u.day).decompose()

# * Manual input
# min_alt = coord.Angle("45"+"d")
# min_dec = coord.Angle("30"+"d")
# window_east = coord.Angle("3"+" hours")
# window_west = coord.Angle("1"+" hours")
# obs_date ="2024-02-08"
# obs_utc = "20:00"
# obs_duration = coord.Angle("1"+" hours")
   
# * Date/time/location objects
obs_loc = coord.EarthLocation(lon=obs_lon, lat=obs_lat)
obs_utc = time.Time(obs_date + " " + obs_utc, scale='utc')
obs = time.Time(obs_utc, location=obs_loc)

# * Local sidereal time
st = obs.sidereal_time("apparent") # ? apparent or absolute

# * Borders and corners coordinates
dec_upper = min_dec
dec_lower = + min_alt + obs_loc.lat - coord.Angle("90d") 
a_west = st + obs_sky_rotation + window_west
a_east = st - window_east

# * Output
print("")
print("\033[36m"+"Location: "+"\033[0m"+ obs_name)
print("\033[36m"+"\tlat: "+"\033[0m" + obs_lat.to_string())
print("\033[36m"+"\tlon: "+"\033[0m" + obs_lon.to_string())
print("\033[36m"+"Observation "+"\033[0m")
print("\033[36m"+"\tdate and time: "+"\033[0m"+ obs_utc.to_string())
print("\033[36m"+"\tsidereal time: "+"\033[0m"+ st.to_string(unit=u.hour))
print("\033[36m"+"\tduration:      "+"\033[0m"+ obs_duration.to_string(unit=u.hour))
print("\033[36m"+"\tsky rotation:  "+"\033[0m"+ obs_sky_rotation.to_string(unit=u.hour))
print("\033[36m"+"Sky coordinates window:")
print("\033[36m"+"\teast ra:   "+"\033[0m" + a_east.to_string(unit=u.hour))
print("\033[36m"+"\twest ra:   "+"\033[0m" + a_west.to_string(unit=u.hour))
print("\033[36m"+"\tupper dec: "+"\033[0m" + dec_upper.to_string(unit=u.degree))
print("\033[36m"+"\tlower dec: "+"\033[0m" + dec_lower.to_string(unit=u.degree))

print("\033[36m"+"Simbad constraints:"+"\033[0m")
print("\033[36m"+"\tRA: "+"\033[0m" + "<" + a_west.to_string(unit=u.hour, sep=":") + " && " + ">" + a_east.to_string(unit=u.hour, sep=":"))
print("\033[36m"+"\tDE: "+"\033[0m" + "<" + dec_upper.to_string(unit=u.degree, sep=":") + " && " + ">" + dec_lower.to_string(unit=u.degree, sep=":"))

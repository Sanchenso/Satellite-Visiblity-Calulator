from astropy.coordinates import EarthLocation, GCRS, ITRS, AltAz, CartesianRepresentation
from astropy.time import Time
import astropy.units as u

# Исходные данные
eci_ka = [4435144, -2137297, 4670064] * u.m   # положение КА в ECI (GCRS)
JDN = 8084.185608609847                        # дни от J2000 (TT)
observer_lat = 45.920266 * u.deg
observer_lon = -63.342286 * u.deg
h = 0 * u.m
min_elevation = 15 * u.deg

# расчет времени
t = Time(2451545.0 + JDN, format='jd', scale='tt')  # JD(TT) от J2000
t_ut1 = t.ut1  # дельта TT-UT1 из astropy

# Координаты наблюдателя
observer = EarthLocation(lat=observer_lat, lon=observer_lon, height=h)

# Перевод положения КА из ECI (GCRS) в топоцентрические координаты AltAz
gcrs = GCRS(CartesianRepresentation(eci_ka), obstime=t)
itrs = gcrs.transform_to(ITRS(obstime=t_ut1))
altaz = itrs.transform_to(AltAz(obstime=t, location=observer))

elevation = altaz.alt
print(f"Угол возвышения: {elevation:.2f}, {'КА виден' if elevation > min_elevation else 'КА не виден'}")
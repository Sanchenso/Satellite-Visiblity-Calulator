"""
Программа для оценки видимости космического аппарата (КА) с точки зрения наземного наблюдателя.
Алгоритм выполняет пошаговое преобразование координат спутника из инерциальной системы координат (ECI)
в земную систему (ECEF) с учетом вращения Земли и прецессии-нутации по модели IAU 2000B
"""


import math
import numpy as np
import erfa

# =============================================================================
# Константы и базовые параметры Земли по WGS-84
# =============================================================================

F = 1/298.257223563    # сплюснутость, данные взяты из "Руководство по Всемирной геодезической системе — 1984 (WGS-84)"
A = 6378137            # Большая полуось эллипсоида, м
E2 = F * (2 - F)       # квадрат эксцентриситета эллипсоида
day_sec = 86400.0      # секунд в сутках


# =============================================================================
# Вспомогательные функции 
# =============================================================================
def jd_from_tt_offset(JDN_tt_from_J2000):
    """
    Преобразование дней ТТ от эпохи J2000 в полный юлианский день
    """
    return 2451545.0 + JDN_tt_from_J2000

def tt_utc_leap_seconds(leap_seconds):
    """
    перевод TT в UTC, сек
    """
    return 32.184 + float(leap_seconds) 

def era_from_ut1(jd_ut1):
    """
    расчет угла вращения Земли (ERA) по UT1, модель IAU/IERS (CIO-подход)
    возвращает угол ERA, рад
    """
    D_ut1 = (jd_ut1 - 2451545.0)                               # полный юлианский день UT1 от 2000-01-01 12:00
    f = (0.7790572732640 + 1.00273781191135448 * D_ut1) % 1.0  # дробная часть оборота
    return 2 * math.pi * f 

def geodetic_to_ecef(lat_deg, lon_deg, h_m):
    """
    # Перевод геодезических координат наблюдателя в геоцентрическую систему (ECEF) по WSG-84
    # возвращает координаты ECEF
    """
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    N = A / math.sqrt(1.0 - E2 * math.sin(lat)**2) # радиус кривизны в главном вертикале
    X = (N + h_m) * math.cos(lat) * math.cos(lon)  # преобразование координат наблюдателя из геодезических в ECEF
    Y = (N + h_m) * math.cos(lat) * math.sin(lon)
    Z = (N * (1.0 - E2) + h_m) * math.sin(lat)
    return (X, Y, Z)

def enu_components(dxyz, lat_deg, lon_deg):
    """
    Преобразование вектора из геоцентрической ECEF в топоцентрическую проекцию East(E), North(N), Up(U)
    """
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sinp, cosp = math.sin(lat), math.cos(lat)
    sinl, cosl = math.sin(lon), math.cos(lon)
    dx, dy, dz = dxyz
    E = -sinl*dx + cosl*dy
    N = -sinp*cosl*dx - sinp*sinl*dy + cosp*dz
    U =  cosp*cosl*dx + cosp*sinl*dy + sinp*dz
    return E, N, U

def elevation_from_enu(E, N, U):
    """
    Возвращает угол возвышения из ENU, в град
    """
    return math.degrees(math.atan2(U, math.hypot(E, N)))

# =============================================================================
# Функции вращения
# =============================================================================

def R1(x):
    """
    вращение вокруг оси X
    """
    cx, sx = math.cos(x), math.sin(x)
    return np.array([[1, 0, 0],
                     [0, cx, -sx],
                     [0, sx,  cx]])

def R2(y):
    """
    вращение вокруг оси Y
    """
    cy, sy = math.cos(y), math.sin(y)
    return np.array([[ cy, 0,  sy],
                     [  0, 1,   0],
                     [-sy, 0,  cy]])

def R3(z):
    """
    вращение вокруг оси Z
    """
    cz, sz = math.cos(z), math.sin(z)
    return np.array([[ cz,  sz, 0],
                     [-sz,  cz, 0],
                     [  0,   0, 1]])

def W_polar(xp, yp):
    """
    построение матрицы полярного движения W(xp, yp)
    параметры смещения оси вращения Земли xp, yp берутся из EOP, см. IERS https://iers-conventions.obspm.fr/content/tn36.pdf?utm_source=chatgpt.com
    """
    return R2(xp) @ R1(yp)


def C_from_XYs(jd_tt):
    """
    Матрица прецессии-нутации CIP/CIO по модели укороченного ряда IAU 2000B (через X, Y, s)
    Построение матрицы C с помощью Erfa, из ECI в CIRS
    """
    jd1, jd2 = jd_tt, 0.0
    X, Y, s = erfa.xys00b(jd1, jd2)  # X, Y координаты полюса; s - угол, определяющий ориентацию оси X CIRS
    C = erfa.c2ixys(X, Y, s)         # матрица C 3x3
    return C

# =============================================================================
# Преобразование ECI в ECEF
# =============================================================================


def eci_to_ecef_gcrs(r_gcrs, jd_tt, leap_seconds=37, DUT1=0.0, xp_as=0.0, yp_as=0.0):
    """
    Преобразование координат спутника из инерциальной системы ECI
    в земную систему координат (ECEF), по стандарту IAU/IERS

    Параметры:
        r_gcrs - координаты спутника в GCRS, м
        jd_tt - юлианская дата в шкале TT
        leap_seconds - количество високосных секунд (TAI−UTC)
        DUT1 - разность UT1-UTC, сек
        xp_as, yp_as - параметры полярного движения, угл. сек

    Возвращает координаты спутника в ECEF, м
    """

    #  перевод TT в UTC
    tt_minus_utc = tt_utc_leap_seconds(leap_seconds) / day_sec
    jd_utc = jd_tt - tt_minus_utc
    #  перевод UTC в UT1
    jd_ut1 = jd_utc + DUT1 / day_sec

    # угол вращения Земли (ERA) по UT1
    theta = era_from_ut1(jd_ut1)

    # Вращение Земли
    R = R3(theta)

    # Прецессия-нутация 
    C = C_from_XYs(jd_tt)  

    # Полярное движение (из IERS Bulletin A)
    xp = xp_as * erfa.DAS2R
    yp = yp_as * erfa.DAS2R
    W = W_polar(xp, yp)

    # Полная матрица GCRS в ECEF
    C2T = W @ R @ C
    return C2T @ r_gcrs

# =============================================================================
# Основной расчет
# =============================================================================

def main():
    # Исходные
    eci_gcrs = np.array([4435144.0, -2137297.0, 4670064.0])  # Координаты КА в ECI
    min_elevation = 15
    lat_deg= 45.920266      # северной широты
    lon_deg = -63.342286    # западной долготы
    h = 0                   # высота над уровнем моря
    JDN = 8084.185608609847 # номер юлианского дня
    leap_seconds = 37
    jd_tt = jd_from_tt_offset(JDN)

    # EOP на дату
    DUT1 = 0.0    # сек
    xp_as = 0.0   # угл сек
    yp_as = 0.0   # угл сек

    # Спутник, перевод из GCRS в ECEF:
    r_ecef = eci_to_ecef_gcrs(eci_gcrs, jd_tt, leap_seconds=leap_seconds, DUT1=DUT1, xp_as=xp_as, yp_as=yp_as)

    # Наблюдатель, перевод из геодезической в ECEF
    r_site = geodetic_to_ecef(lat_deg, lon_deg, h)

    # Топоцентрическая проекция ENU 
    dxyz = r_ecef - r_site
    E, N, U = enu_components(dxyz, lat_deg, lon_deg)

    # Угол возвышения
    elevation = elevation_from_enu(E, N, U)
    visible = min_elevation <= elevation

    # вывод в консоль
    print(f"Угол возвышения: {elevation:.2f}°, {'КА виден' if visible else 'КА не виден'}")

if __name__ == "__main__":
    main()
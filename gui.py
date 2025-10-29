import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import sat_math

def parse_float(entry, name, lo=None, hi=None):
    """Проверка ввода"""
    s = entry.get().strip().replace(",", ".")
    try:
        x = float(s)
    except ValueError:
        raise ValueError(f"Поле «{name}» должно быть числом.")
    if lo is not None and x < lo:
        raise ValueError(f"Поле «{name}» должно быть ≥ {lo}.")
    if hi is not None and x > hi:
        raise ValueError(f"Поле «{name}» должно быть ≤ {hi}.")
    return x

def parse_int(entry, name, lo=None, hi=None):
    """Проверка ввода"""
    s = entry.get().strip()
    try:
        v = int(s)
    except ValueError:
        raise ValueError(f"Поле «{name}» должно быть целым числом.")
    if lo is not None and v < lo:
        raise ValueError(f"Поле «{name}» должно быть ≥ {lo}.")
    if hi is not None and v > hi:
        raise ValueError(f"Поле «{name}» должно быть ≤ {hi}.")
    return v

def set_eop_fields_state():
    state = "normal" if var_use_eop.get() else "disabled"
    for w in (ent_dut1, ent_xp, ent_yp):
        w.configure(state=state)

def compute():
    try:
        x = parse_float(ent_x, "ECI X (м)")
        y = parse_float(ent_y, "ECI Y (м)")
        z = parse_float(ent_z, "ECI Z (м)")
        jdn_tt = parse_float(ent_jdn, "JDN от J2000 (TT)")

        lat = parse_float(ent_lat, "Широта (°)", lo=-90.0, hi=90.0)
        lon = parse_float(ent_lon, "Долгота (°)", lo=-180.0, hi=180.0)
        h   = parse_float(ent_h, "Высота (м)")

        min_elev = parse_float(ent_min_elev, "Порог видимости (°)", lo=-90.0, hi=90.0)

        leap_seconds = parse_int(ent_leaps, "Leap seconds", lo=0, hi=100)

        if var_use_eop.get():
            dut1 = parse_float(ent_dut1, "DUT1 (с)", lo=-1.0, hi=1.0)  # разумные рамки
            xp_as = parse_float(ent_xp, "xp (″)", lo=-1.0, hi=1.0)
            yp_as = parse_float(ent_yp, "yp (″)", lo=-1.0, hi=1.0)
        else:
            dut1 = 0.0
            xp_as = 0.0
            yp_as = 0.0

        r_gcrs = np.array([x, y, z], dtype=float)
        jd_tt = 2451545.0 + jdn_tt  # полный JD(TT)

        r_itrs = sat_math.eci_to_ecef_gcrs(
            r_gcrs, jd_tt,
            leap_seconds=leap_seconds,
            DUT1=dut1,
            xp_as=xp_as,
            yp_as=yp_as
        )

        r_site = sat_math.geodetic_to_ecef(lat, lon, h)

        dxyz = r_itrs - r_site
        E, N, U = sat_math.enu_components(dxyz, lat, lon)
        elev = sat_math.elevation_from_enu(E, N, U)
        visible = (elev >= min_elev)

        lbl_result.config(
            text=(
                f"КА ECEF (м): X={r_itrs[0]:.2f}, Y={r_itrs[1]:.2f}, Z={r_itrs[2]:.2f}\n"
                f"Угол возвышения: {elev:.2f}° — {'КА виден' if visible else 'КА не виден'}"
            )
        )
    except Exception as e:
        messagebox.showerror("Ошибка ввода", str(e))

#GUI
root = tk.Tk()
root.title("Visibility (ECI→ECEF) — demo")

frm = ttk.Frame(root, padding=12)
frm.grid(sticky="nsew")
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

def add_row(r, label, default="", width=16):
    ttk.Label(frm, text=label).grid(row=r, column=0, sticky="w", padx=(0,6), pady=3)
    e = ttk.Entry(frm, width=width)
    e.insert(0, default)
    e.grid(row=r, column=1, sticky="ew", pady=3)
    return e

# Вход в ECI
ent_x = add_row(0, "ECI X (м):", "4435144")
ent_y = add_row(1, "ECI Y (м):", "-2137297")
ent_z = add_row(2, "ECI Z (м):", "4670064")

ttk.Separator(frm).grid(row=3, column=0, columnspan=2, sticky="ew", pady=4)

# Входные данные по умолчанию
ent_jdn = add_row(4, "JDN от J2000 (TT):", "8084.185608609847")
ent_lat = add_row(5, "Широта (°):", "45.920266")
ent_lon = add_row(6, "Долгота (°):", "-63.342286")
ent_h   = add_row(7, "Высота (м):", "0")
ent_min_elev = add_row(8, "Порог видимости (°):", "15")

ttk.Separator(frm).grid(row=9, column=0, columnspan=2, sticky="ew", pady=4)

# Время/ориентация Земли
ent_leaps = add_row(10, "Leap seconds:", "37")

# Переключатель EOP
var_use_eop = tk.BooleanVar(value=False)
chk = ttk.Checkbutton(frm, text="Учитывать EOP (DUT1, xp, yp)", variable=var_use_eop, command=set_eop_fields_state)
chk.grid(row=11, column=0, columnspan=2, sticky="w", pady=(2,2))

ent_dut1  = add_row(12, "DUT1 (с):", "0")
ent_xp    = add_row(13, "xp (″):", "0")
ent_yp    = add_row(14, "yp (″):", "0")

set_eop_fields_state()

# Кнопка расчёта
btn = ttk.Button(frm, text="Рассчитать", command=compute)
btn.grid(row=15, column=0, columnspan=2, sticky="ew", pady=(8,6))

# Результат
lbl_result = ttk.Label(frm, text="—", justify="left")
lbl_result.grid(row=16, column=0, columnspan=2, sticky="w", pady=(6,0))

for c in (0,1):
    frm.columnconfigure(c, weight=1)

root.mainloop()
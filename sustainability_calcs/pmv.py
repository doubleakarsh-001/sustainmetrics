import math

def calculate_pmv(temperature_c, humidity, wind_speed, met=1.2, clo=0.5, tr=None, pa=None):
    """
    Calculate the Predicted Mean Vote (PMV) based on ISO 7730.
    Args:
        temperature_c (float): Air temperature in Celsius.
        humidity (float): Relative humidity in percent (0-100).
        wind_speed (float): Air velocity in m/s.
        met (float): Metabolic rate (default 1.2 met).
        clo (float): Clothing insulation (default 0.5 clo).
        tr (float): Mean radiant temperature in Celsius (default: same as air temperature).
        pa (float): Partial water vapor pressure in kPa (default: calculated from humidity).
    Returns:
        float: PMV value.
    """
    # Constants
    M = met * 58.15  # metabolic rate in W/m^2
    W = 0            # external work in W/m^2
    Icl = clo * 0.155  # clothing insulation in m^2K/W
    ta = temperature_c
    if tr is None:
        tr = ta
    # Calculate partial water vapor pressure (pa) if not provided
    if pa is None:
        # Saturation vapor pressure (kPa)
        pws = 0.6108 * math.exp(17.27 * ta / (ta + 237.3))
        pa = humidity / 100.0 * pws
    # Heat transfer coefficient by convection
    hc = 2.38 * abs(ta - tr) ** 0.25 if abs(ta - tr) > 0.1 else 12.1 * math.sqrt(wind_speed)
    # Clothing area factor
    fcl = 1.05 + 0.1 * clo if clo > 0.5 else 1 + 0.2 * clo
    # Calculate PMV
    tcl = ta + (35.5 - ta) / (3.5 * (6.45 * Icl + 0.1))
    hl1 = 3.05 * 0.001 * (5733 - (6.99 * (M - W)) - pa * 1000)
    hl2 = 0.42 * ((M - W) - 58.15)
    hl3 = 1.7 * 0.00001 * M * (5867 - pa * 1000)
    hl4 = 0.0014 * M * (34 - ta)
    hl5 = fcl * hc * (tcl - ta)
    hl6 = fcl * 3.96e-8 * ((tcl + 273) ** 4 - (tr + 273) ** 4)
    pmv = (0.303 * math.exp(-0.036 * M) + 0.028) * ((M - W) - hl1 - hl2 - hl3 - hl4 - hl5 - hl6)
    return pmv

# Example usage
if __name__ == "__main__":
    temp = 25  # Celsius
    rh = 50    # %
    wind = 0.1 # m/s
    pmv = calculate_pmv(temp, rh, wind)
    print(f"PMV: {pmv:.2f}")

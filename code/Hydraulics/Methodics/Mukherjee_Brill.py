from Hydraulics.Formulas import *
import warnings
warnings.filterwarnings("ignore")

def calculate (mishenko, tubing):
    """
    :param oil_params: Данные о свойствах нефти и режиме работы скважины
    :param mishenko: Уточненные характеристики нефтегазовой смеси, рассчитанные по методике Мищенко
    :param tubing: Сведенные данные инклинометрии и труб НКТ
    :return: dP/dL - локальная производная давления по длине
    """
    angle = tubing["angle"]



    wm = mishenko.Q * 4 / (math.pi * tubing["IntDiameter"] ** 2)
    g = 9.81
    if (mishenko.VolumeWater == 1):
        Re = wm * tubing["IntDiameter"] * mishenko.CurrentLiquidDensity / mishenko.CurrentWaterViscosity
        lambda_ = 64 / Re if Re < 2300 else 0.11 * (68 / Re + tubing["Roughness"] / tubing["IntDiameter"]) ** 0.25
        return lambda_ * mishenko.CurrentLiquidDensity * (wm ** 2) * 0.5 / tubing["IntDiameter"],  mishenko.CurrentLiquidDensity * mishenko.g

    # Коэффициент по скорости жидкости
    TensionLiquidGas = mishenko.TensionWaterGas * mishenko.VolumeWater + mishenko.TensionOilGas * (1 + mishenko.VolumeWater)
    Lw = mishenko.OilVolumeCoeff * wm * (mishenko.SaturatedOilDensity/ (g * TensionLiquidGas)) ** 0.25
    # Коэффициент по скорости свободного газа
    Lg = mishenko.VolumeGas * wm * (mishenko.SaturatedOilDensity / (g * TensionLiquidGas)) ** 0.25
    # Коэффициент по вязкости жидкости
    Lm = mishenko.CurrentOilViscosity * wm * (mishenko.SaturatedOilDensity / (g * TensionLiquidGas ** 3)) ** 0.25
    mu = mishenko.OilVolumeCoeff * mishenko.CurrentOilViscosity + mishenko.VolumeGas * mishenko.CurrentFreeGasViscosity
    dens = mishenko.CurrentLiquidDensity + mishenko.VolumeGas * mishenko.CurrentFreeGasDensity
    # Число Рейнольдса
    Re = dens * wm * tubing["IntDiameter"] / mu
    if (mishenko.CurrentP >= mishenko.SaturationPressure):
        Re = wm * tubing["IntDiameter"] * mishenko.CurrentLiquidDensity / mishenko.CurrentOilViscosity
        lambda_ = 64 / Re if Re < 2300 else 0.11 * (68 / Re + tubing["Roughness"] / tubing["IntDiameter"]) ** 0.25
        return lambda_ * mishenko.CurrentLiquidDensity * wm ** 2 * 0.5 / tubing["IntDiameter"], mishenko.CurrentLiquidDensity * mishenko.g
    try:
        form = get_flow_structure_MB(Lw, Lg, Lm, tubing)
    except:
        form = "stratified"
    C1, C2, C3, C4, C5, C6 = get_coefs_MB(tubing, form)
    #!!!!!!!!!!!
    try:
        phi1 = math.exp((C1 + C2 * math.sin(math.radians(angle)) + C3 * math.cos(
        math.radians(angle)) ** 2 + C4 * Lm ** 2) * (Lg ** C5 / Lw ** C6))
    except:
        phi1 = 1
    phi1 = min(phi1, 1)
    # Число Фруда смеси
    Fr = count_Frud(mishenko, wm, tubing)
    # Объемные концентрации
    wg = 3.3 * (g * mishenko.TensionOilGas / (mishenko.SaturatedOilDensity - mishenko.CurrentFreeGasDensity)) ** 0.25 * (
            mishenko.SaturatedOilDensity / mishenko.CurrentFreeGasDensity) ** 0.5
    # Относительная скорость реверса
    u = wm / wg
    a, b, kr, F = get_coefficients(u,  mishenko.CurrentOilViscosity, mishenko.CurrentFreeGasViscosity, mishenko.SaturatedOilDensity,
                                   mishenko.CurrentFreeGasDensity, Fr)
    # Истинная объемная концентрация
    phi2 = mishenko.VolumeGas * (1 - (1 - kr) * (a - mishenko.VolumeGas) / (b - mishenko.VolumeGas))
    dens_true = mishenko.CurrentLiquidDensity
    # Коэффициент расширения потока
    Ek = mishenko.VolumeGas * wm ** 2 * dens / mishenko.CurrentP if mishenko.CurrentP >=1 else mishenko.VolumeGas * wm ** 2 * dens
    # Коэффициент гидравлического сопротивления однофазного потока
    if Re < 2000:
        lambda0 = 64 / Re
    else:
        lambda0 = (-2 * math.log10(2 * tubing["Roughness"] / (3.7 * tubing["IntDiameter"]) - 5.02 * math.log10(
            2 * tubing["Roughness"] / (3.7 * tubing["IntDiameter"]) + 13 / Re))) ** -2
    #Локальный градиент давления
    dP_fric, dP_grav = count_dP_MB(mishenko, tubing, form, lambda0, dens_true, Ek, phi1, phi2, Lw, Lg, Lm, wm)

    return dP_fric, dP_grav

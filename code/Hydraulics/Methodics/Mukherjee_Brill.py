from Hydraulics.Formulas import *
import warnings
from Hydraulics.Properties.Mishenko import Mishenko
warnings.filterwarnings("ignore")

def calc_lambda(Re, tubing):
    if Re < 2000:
        lambda0 = 64 / Re
    else:
        lambda0 = 0.11 * (tubing["Roughness"] / tubing["IntDiameter"] + 68.5 / Re) ** 0.25

    for i in range(10):
        lambda0 = (1.74 - 2 * math.log10(2 * tubing["Roughness"] / (tubing["IntDiameter"] ** 2) + 18.7 / (Re * (lambda0 ** 0.5)))) ** -2
    return lambda0


def calculate (mishenko : Mishenko, tubing):
    """
    :param oil_params: Данные о свойствах нефти и режиме работы скважины
    :param mishenko: Уточненные характеристики нефтегазовой смеси, рассчитанные по методике Мищенко
    :param tubing: Сведенные данные инклинометрии и труб НКТ
    :return: dP/dL - локальная производная давления по длине
    """
    angle = tubing["angle"]

    wm = mishenko.Q_liq_and_gas_m3_s * 4 / (math.pi * tubing["IntDiameter"] ** 2)
    g = 9.81
    if (mishenko.VolumeWater_fraction == 1):
        Re = wm * tubing["IntDiameter"] * mishenko.CurrentLiquidDensity_kg_m3 / mishenko.CurrentWaterViscosity_Pa_s
        lambda_ = calc_lambda(Re, tubing)
        return lambda_ * mishenko.CurrentLiquidDensity_kg_m3* (wm ** 2) * 0.5 / tubing["IntDiameter"],  mishenko.CurrentLiquidDensity_kg_m3 * mishenko.g

    # Коэффициент по скорости жидкости
    TensionLiquidGas = mishenko.TensionWaterGas * mishenko.VolumeWater_fraction + mishenko.TensionOilGas * (1 + mishenko.VolumeWater_fraction)
    TensionLiquidGas = abs(TensionLiquidGas)
    Lw = mishenko.oil_params.volumeoilcoeff * wm * (mishenko.SaturatedOilDensity_kg_m3/ (g * TensionLiquidGas)) ** 0.25
    # Коэффициент по скорости свободного газа
    Lg = mishenko.VolumeGas_fraction * wm * (mishenko.SaturatedOilDensity_kg_m3 / (g * TensionLiquidGas)) ** 0.25
    # Коэффициент по вязкости жидкости
    Lm = mishenko.CurrentOilViscosity_Pa_s * wm * (mishenko.SaturatedOilDensity_kg_m3 / (g * TensionLiquidGas ** 3)) ** 0.25
    mu = mishenko.oil_params.volumeoilcoeff * mishenko.CurrentOilViscosity_Pa_s + mishenko.VolumeGas_fraction * mishenko.CurrentFreeGasViscosity_Pa_s
    dens = mishenko.CurrentLiquidDensity_kg_m3 * (1 - mishenko.VolumeGas_fraction) + mishenko.VolumeGas_fraction * mishenko.CurrentFreeGasDensity_kg_m3
    # Число Рейнольдса
    Re = dens * wm * tubing["IntDiameter"] / mu
    if (mishenko.CurrentP_MPa >= mishenko.SaturationPressure_MPa):
        Re = wm * tubing["IntDiameter"] * mishenko.CurrentLiquidDensity_kg_m3 / mishenko.CurrentOilViscosity_Pa_s
        lambda_ = calc_lambda(Re, tubing)
        return lambda_ * mishenko.CurrentLiquidDensity_kg_m3 * wm ** 2 * 0.5 / tubing["IntDiameter"], mishenko.CurrentLiquidDensity_kg_m3 * mishenko.g
    try:
        form = get_flow_structure_MB(Lw, Lg, Lm, tubing)
    except:
        form = "stratified"
    C1, C2, C3, C4, C5, C6 = get_coefs_MB(tubing, form)
    #!!!!!!!!!!!
    try:
        phi1 = math.exp((C1 + C2 * math.sin(math.radians(angle)) + C3 * math.cos(math.radians(angle)) ** 2 + C4 * Lm ** 2) * (Lg ** C5 / Lw ** C6))
    except:
        phi1 = 1
    phi1 = min(phi1, 1)
    # Число Фруда смеси
    Fr = count_Frud(mishenko, wm, tubing)
    # Объемные концентрации
    wg = 3.3 *abs (g * mishenko.TensionOilGas / (mishenko.SaturatedOilDensity_kg_m3 - mishenko.CurrentFreeGasDensity_kg_m3)) ** 0.25 * abs(
            mishenko.SaturatedOilDensity_kg_m3 / mishenko.CurrentFreeGasDensity_kg_m3) ** 0.5
    # Относительная скорость реверса
    u = wm / wg
    a, b, kr, F = get_coefficients(u,  mishenko.CurrentOilViscosity_Pa_s, mishenko.CurrentFreeGasViscosity_Pa_s, mishenko.SaturatedOilDensity_kg_m3,
                                   mishenko.CurrentFreeGasDensity_kg_m3, Fr)
    # Истинная объемная концентрация
    phi2 = mishenko.VolumeGas_fraction * (1 - (1 - kr) * (a - mishenko.VolumeGas_fraction) / (b - mishenko.VolumeGas_fraction))
    dens_true = mishenko.CurrentLiquidDensity_kg_m3
    # Коэффициент расширения потока
    Ek = mishenko.VolumeGas_fraction * wm ** 2 * mishenko.CurrentFreeGasDensity_kg_m3 / (mishenko.CurrentP_MPa * 1e6) #if mishenko.CurrentP_MPa >=1 else mishenko.VolumeGas_fraction * wm ** 2 * mishenko.CurrentFreeGasDensity_kg_m3
    Ek = min(Ek, 0.5)
    # Коэффициент гидравлического сопротивления однофазного потока

    lambda0 = calc_lambda(Re, tubing)

    #Локальный градиент давления
    dP_fric, dP_grav = count_dP_MB(mishenko, tubing, form, lambda0, dens_true, Ek, phi1, phi2, Lw, Lg, Lm, wm)

    return dP_fric, dP_grav

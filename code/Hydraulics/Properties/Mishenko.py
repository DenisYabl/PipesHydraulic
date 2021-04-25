import pandas as pd
import math
from Hydraulics.Formulas import get_dens_freegas
from Tools.HE2_ABC import oil_params
from Tools.HE2_Logger import check_for_nan, getLogger
from collections import namedtuple

# TODO: Диаметр трубы нужен для уточнения вязкости смеси при движении по трубе, больше ни для чего
# Поэтому надо разделить расчет на две части - когда поток идет по трубе, и когда нет (например через насос)
# В первом случае мы считаем вязкость через режим течения и проскальзывание
# Во втором либо линейной интерполяцией, либо инверсией фаз

field_list = [
'oil_params',
'VolumeWater_fraction',
"Q_liq_m3_s",
"g",
'SaturationPressure_MPa',
'CurrentP_MPa',
'CurrentT_K',
"DissolvedGasAmount",
'FreeGasDensity_kg_m3',
'SaturatedOilDensity_kg_m3',
'CurrentWaterDensity_kg_m3',
'CurrentOilViscosity_Pa_s',
'CurrentWaterViscosity_Pa_s',
'RelativePressure',
'RelativeTemperature',
'CurrentFreeGasDensity_kg_m3',
'CurrentFreeGasViscosity_Pa_s',
'CurrentLiquidDensity_kg_m3',
'TensionOilGas',
'TensionWaterGas',
'TensionOilWater',
'Q_gas_m3_s',
'Q_liq_and_gas_m3_s',
'VolumeGas_fraction',
'Thermal_capacity'

]

Mishenko = namedtuple('Mishenko', field_list)

# def print(mishenko : Mishenko):
#     print(f"""Газовый фактор GasFactor= {mishenko.oil_params.gasFactor}
# Обводненность VolumeWater =  {mishenko.oil_params.VolumeWater}
# Дебит mishenko.Q =  {mishenko.Q}
# Объемный фактор нефти OilVolumeCoeff {mishenko.OilVolumeCoeff}
# Текущее давление насыщения SaturationPressure_MPa = {mishenko.SaturationPressure_MPa}
# Текущее давление CurrentP = {mishenko.CurrentP}
# Текущая температуры CurrentT = {mishenko.CurrentT}
# Удельный объем растворенного газа DissolvedGasAmount = {mishenko.DissolvedGasAmount}
# Плотность газонасыщенной нефти SaturatedOilDensity = {mishenko.SaturatedOilDensity}
# Плотность воды в текущих условиях CurrentWaterDensity = {mishenko.CurrentWaterDensity}
# Вязкость нефти в текущих условиях CurrentOilViscosity = {mishenko.CurrentOilViscosity}
# Вязкость воды в текущих условиях CurrentWaterViscosity = {mishenko.CurrentWaterViscosity}
# Плотность свободного газа в текущих условиях CurrentFreeGasDensity = {mishenko.CurrentFreeGasDensity}
# Вязкость газа в текущих условиях CurrentFreeGasViscosity = {mishenko.CurrentFreeGasViscosity}
# Плотность жидкой фазы в текущих условиях CurrentLiquidDensity = {mishenko.CurrentLiquidDensity}
# Объемная доля газа в смеси VolumeGas = {mishenko.VolumeGas}""")



def from_oil_params(P_bar, T_C, X_kg_sec, calc_params:oil_params, tubing=None):
    SaturationPressure_MPa = calc_params.sat_P_bar * 101325 * 1e-6
    CurrentP = P_bar * 101325 * 1e-6
    PlastT = calc_params.plastT_C + 273
    CurrentT = T_C + 273
    GasFactor = calc_params.gasFactor
    Saturation_pressure = SaturationPressure_MPa - (PlastT - CurrentT) / (GasFactor * (0.91 - 0.09))
    if (CurrentP >= Saturation_pressure) | (calc_params.volumewater_percent == 100):
        return two_phase_flow(P_bar, T_C, X_kg_sec, calc_params, tubing)
    else:
        return three_phase_flow(P_bar, T_C, X_kg_sec, calc_params)

def two_phase_flow(P_bar, T_C, X_kg_sec, calc_params:oil_params, tubing=None):
    """
    :param oil_params: Параметры нефти
    """
    # Давление насыщения нефти попутным газом
    SaturationPressure_MPa = calc_params.sat_P_bar * 101325 * 1e-6
    # Текущее давление
    CurrentP = P_bar * 101325 * 1e-6
    # Пластовая температура
    PlastT = calc_params.plastT_C + 273
    # Температура в текущей точке
    CurrentT = T_C + 273
    # Газовый фактор нефти
    GasFactor = calc_params.gasFactor
    # Доля углеводородных газов
    CarbGasAmount = 0.91  # Нет в данных
    # Доля неуглеводородных газов
    NonCarbGasAmount = 0.09  # Нет в данных
    # Плотность дегазированной нефти
    SepOilDensity = calc_params.oildensity_kg_m3
    # Плотность газа
    GasDensity = calc_params.gasdensity_kg_m3
    # Динамическая вязкость сепарированной нефти
    SepOilDynamicViscosity = calc_params.oilviscosity_Pa_s
    # Обводненность нефти
    VolumeWater = calc_params.volumewater_percent / 100
    # Плотность пластовой воды
    PlastWaterDensity = calc_params.waterdensity_kg_m3
    # Объемный фактор нефти
    OilVolumeCoeff = calc_params.volumeoilcoeff

    g = 9.81

    Saturation_pressure = SaturationPressure_MPa - (PlastT - CurrentT) / (GasFactor * (0.91 - 0.09))

    #Объемное расходное водосодержание
    VolumeWater = VolumeWater / (VolumeWater + OilVolumeCoeff * (1 - VolumeWater))

    #Поверхностное натяжение нефть-газ
    TensionOilGas = 10 ** -(1.58 + 0.05 * CurrentP) - 72e-6 * (CurrentT - 303)

    # Поверхностное натяжение вода - газ
    TensionWaterGas = 10 ** -(1.19 + 0.01 * CurrentP)

    # Поверхностное натяжение нефть - вода
    TensionOilWater = TensionWaterGas - TensionOilGas

    # Вязкость нефти в текущих условиях
    if SepOilDensity > 860:
        an = 1e-3 * (2.513 - 0.001975 * SepOilDensity)
    else:
        an = 1e-3 * (3.083 - 0.002638 * SepOilDensity)

    vzgr = 0.001055 * (1 + 5 * an) * SepOilDynamicViscosity * SepOilDensity

    C = 1 / (1 + 0.00144 * (CurrentT - 293) * math.log(1000 * SepOilDynamicViscosity))
    Amnt = 1e-3 * (1000 * SepOilDynamicViscosity) ** C

    A = 1 + 0.0129 * vzgr - 0.0364 * vzgr ** 0.85
    B = 1 + 0.0017 * vzgr - 0.0228 * vzgr ** 0.667

    CurrentOilViscosity = A * Amnt ** B

    # Вязкость воды в текущих условиях

    CurrentWaterViscosity = 0.00178 / (1 + 0.037 * (CurrentT - 273) + 0.00022 * (CurrentT - 273) ** 2)

    # Вспомогательные коэффициенты
    a1 = 1 + 0.0054 * (CurrentT - 303)
    a2 = (1 + math.log(abs(CurrentP), 10) / (1 + math.log(abs(SaturationPressure_MPa), 10)) - 1)  # !!!!!!!!!!!
    a3 = SepOilDensity * GasFactor * 1e-3 - 186
    a4 = (3.083 - 2.638e-3 * GasDensity) * 1e-3

    FreeGasDensity = get_dens_freegas(GasDensity, a1, a2, a3)
    # Плотность газонасыщенной нефти
    m = 1 + 0.029 * (CurrentT - 303) * (SepOilDensity * FreeGasDensity * 1e-3 - 0.7966)
    # Относительная плотность растворенного в нефти свободного газа
    if CurrentP < 0:
        DissolvedGasDensity = 0
    else:
        DissolvedGasDensity = GasFactor * (a1 * m * GasDensity) / GasFactor

    SaturatedOilDensity = SepOilDensity * (1 + 1.293e-3 * DissolvedGasDensity * GasFactor /
                                             (a1 * m)) / OilVolumeCoeff
    VolumeOil = 1 - VolumeWater
    #CurrentLiquidDensity = VolumeWater * PlastWaterDensity + VolumeOil * SaturatedOilDensity
    CurrentLiquidDensity = VolumeWater * PlastWaterDensity + VolumeOil * SaturatedOilDensity

    if not X_kg_sec:
        X_kg_sec = 0
    # Объемный расход нефтеводогазовой смеси в условиях транспорта
    Q = X_kg_sec / CurrentLiquidDensity

    #Скорость смеси
    structure = 'emulsion'
    if tubing:
        wm = Q * 4 / (math.pi * tubing["IntDiameter"] ** 2)
        #Критическая скорость смеси
        wkr = 0.457 * math.sqrt(g * tubing["IntDiameter"])
        if wm < wkr:
            structure = "drop"

    if structure == "drop":
        mixture_type = "oil_water" if VolumeWater > 0.5 else "water_oil"
        if mixture_type == "water_oil":
            CurrentOilViscosity = CurrentOilViscosity
        else:
            CurrentOilViscosity = CurrentWaterViscosity

    elif structure == 'emulsion' and tubing:
        wm = Q * 4 / (math.pi * tubing["IntDiameter"] ** 2)
        shiftvelocity = 8 * wm / tubing["IntDiameter"]
        A = (1 + 20 * VolumeWater ** 2) / shiftvelocity ** (0.48 * VolumeWater)
        B = CurrentWaterViscosity if A <= 1 else A * CurrentWaterViscosity
        if VolumeWater > 0.5:
            CurrentOilViscosity = B * (1 + 2.9 * VolumeWater) / (1 - VolumeWater)
    else:
        assert tubing is None

    CurrentOilViscosity = CurrentWaterViscosity if VolumeWater == 1 else CurrentOilViscosity

    Thermal_capacity = VolumeWater * 4183 + (1 - VolumeWater) * 2100
    return Mishenko(oil_params=calc_params, CurrentP_MPa=CurrentP, CurrentT_K=CurrentT, VolumeWater_fraction=VolumeWater, Q_liq_m3_s=Q,
                    g=g, SaturationPressure_MPa=Saturation_pressure,
                    DissolvedGasAmount=1, FreeGasDensity_kg_m3=0,
                    SaturatedOilDensity_kg_m3=SaturatedOilDensity,
                    CurrentWaterDensity_kg_m3=PlastWaterDensity, CurrentOilViscosity_Pa_s=CurrentOilViscosity,
                    CurrentWaterViscosity_Pa_s=CurrentWaterViscosity,
                    RelativePressure=-100500, RelativeTemperature=-100500,
                    CurrentFreeGasDensity_kg_m3=0,
                    CurrentFreeGasViscosity_Pa_s=0, CurrentLiquidDensity_kg_m3=CurrentLiquidDensity,
                    TensionOilGas=TensionOilGas,
                    TensionWaterGas=TensionWaterGas, TensionOilWater=TensionOilWater, Q_gas_m3_s=0, Q_liq_and_gas_m3_s=Q,
                    VolumeGas_fraction=0, Thermal_capacity = Thermal_capacity)

def three_phase_flow(P_bar, T_C, X_kg_sec, calc_params):
    """
    :param oil_params: Параметры нефти
    """
    # Давление насыщения нефти попутным газом
    SaturationPressure_MPa = calc_params.sat_P_bar * 101325 * 1e-6
    # Текущее давление
    CurrentP = P_bar * 101325 * 1e-6
    # Пластовая температура
    PlastT = calc_params.plastT_C + 273
    # Температура в текущей точке
    CurrentT = T_C + 273
    # Газовый фактор нефти
    GasFactor = calc_params.gasFactor
    # Доля углеводородных газов
    CarbGasAmount = 0.91  # Нет в данных
    # Доля неуглеводородных газов
    NonCarbGasAmount = 0.09  # Нет в данных
    # Плотность дегазированной нефти
    SepOilDensity = calc_params.oildensity_kg_m3
    # Плотность газа
    GasDensity = calc_params.gasdensity_kg_m3
    # Динамическая вязкость сепарированной нефти
    SepOilDynamicViscosity = calc_params.oilviscosity_Pa_s
    # Обводненность нефти
    VolumeWater = calc_params.volumewater_percent / 100
    # Плотность пластовой воды
    PlastWaterDensity = calc_params.waterdensity_kg_m3
    # Текущийй расход/дебит
    # Объемный фактор нефти
    OilVolumeCoeff = calc_params.volumeoilcoeff

    g = 9.81
    we_know_gas_amount = (pd.notnull(CarbGasAmount) & pd.notnull(NonCarbGasAmount))

    # Давление насыщения при данной температуре
    if we_know_gas_amount:
        Saturation_pressure = SaturationPressure_MPa - (PlastT - CurrentT) / (
                GasFactor * (CarbGasAmount - NonCarbGasAmount))
    else:
        Saturation_pressure = SaturationPressure_MPa - (PlastT - CurrentT) / (GasFactor * (0.91 - 0.09))
    # Количество растворенного в нефти газа
    # Вспомогательный коэффициент
    if we_know_gas_amount:
        t = 0.32 + 1 / (CarbGasAmount * NonCarbGasAmount + 1.567)
    else:
        t = 0.32 + 1 / 1.567
    DissolvedGasAmount = GasFactor * (CurrentP / Saturation_pressure) ** t
    if CurrentP < 0:
        DissolvedGasAmount = 0
    # Количество выделившегося из одного кубометра нефти свободного газа
    FreeGasAmount = GasFactor - DissolvedGasAmount

    # Относительная плотность выделившегося из нефти свободного газа
    # Вспомогательные коэффициенты
    a1 = 1 + 0.0054 * (CurrentT - 303)
    a2 = (1 + math.log(abs(CurrentP), 10) / (1 + math.log(abs(SaturationPressure_MPa), 10)) - 1)  # !!!!!!!!!!!
    a3 = SepOilDensity * GasFactor * 1e-3 - 186
    a4 = (3.083 - 2.638e-3 * GasDensity) * 1e-3

    FreeGasDensity = get_dens_freegas(GasDensity, a1, a2, a3)
    # Плотность газонасыщенной нефти
    m = 1 + 0.029 * (CurrentT - 303) * (SepOilDensity * FreeGasDensity * 1e-3 - 0.7966)
    # Относительная плотность растворенного в нефти свободного газа
    if CurrentP < 0:
        DissolvedGasDensity = 0
    else:
        DissolvedGasDensity = GasFactor * (a1 * m * GasDensity - FreeGasDensity * (
            FreeGasAmount / GasFactor)) / DissolvedGasAmount

    SaturatedOilDensity = SepOilDensity * (1 + 1.293e-3 * DissolvedGasDensity * DissolvedGasAmount /
                                             (a1 * m)) / OilVolumeCoeff
    # Плотность пластовой воды в текущих условиях
    CurrentWaterDensity = PlastWaterDensity / (1 + 0.000444 * (CurrentT - 293))

    # Вязкость нефти в текущих условиях
    if SepOilDensity > 860:
        an = 1e-3 * (2.513 - 0.001975 * SepOilDensity)
    else:
        an = 1e-3 * (3.083 - 0.002638 * SepOilDensity)

    vzgr = 0.001055 * (1 + 5 * an) * SepOilDynamicViscosity * SepOilDensity

    C = 1 / (1 + 0.00144 * (CurrentT - 293) * math.log(1000 * SepOilDynamicViscosity))
    Amnt = 1e-3 * (1000 * SepOilDynamicViscosity) ** C

    A = 1 + 0.0129 * vzgr - 0.0364 * vzgr ** 0.85
    B = 1 + 0.0017 * vzgr - 0.0228 * vzgr ** 0.667

    CurrentOilViscosity = A * Amnt ** B

    # Вязкость воды в текущих условиях

    CurrentWaterViscosity = 0.00178 / (1 + 0.037 * (CurrentT - 273) + 0.00022 * (CurrentT - 273) ** 2)

    # Плотность выделившегося из нефти свободного газа в рабочих условиях

    if we_know_gas_amount:
        densgc = (FreeGasDensity - 0.97 * NonCarbGasAmount) / CarbGasAmount
    else:
        densgc = (FreeGasDensity - 0.97 * 0.91) / 0.09
    P0 = CurrentP / (4.69 - 0.206 * densgc ** 2)
    RelativePressure = P0
    T0 = CurrentT / (97 + 172 * densgc ** 2)
    RelativeTemperature = T0
    T0 = T0 if T0 >= 1.05 else 1.05
    # Коэффициент сверхсжимаемости попутного газа
    if (0 <= P0 <= 3.8) & (1.17 <= T0 <= 2):
        zc = 1 - P0 * (0.18 / (T0 - 0.73) - 0.135) + 0.016 * (P0 ** 3.45) / (T0 ** 6.1)
    elif (0 <= P0 <= 1.45) & (1.05 <= T0 <= 1.17):
        zc = 0.13 * P0 + (6.05 * T0 - 6.26) * T0 / P0 ** 2
    elif (1.45 <= P0 <= 4) & (1.05 <= T0 <= 1.17):
        zc = 1 - P0 * (0.23 + (1.88 - 1.6 * T0) * P0)
    else:
        zc = 1

    if we_know_gas_amount:
        z = zc * CarbGasAmount + NonCarbGasAmount
    else:
        z = zc * 0.91 + 0.09

    CurrentFreeGasDensity = 3521.7 * FreeGasDensity * CurrentP / (z * CurrentT)

    # Вязкость газа в рабочих условиях
    Mg = 28.97 * CurrentFreeGasDensity
    r = P0 / (0.29 * T0)
    E = (CurrentT / T0) ** (1 / 6) * (P0 / CurrentP) ** (2 / 3) * Mg ** -0.5
    mu0 = 0.0101 * (CurrentT - 273) ** 1.8 - 1.07e-1 * math.sqrt(abs(Mg))  # !!!!!!!!!!!!!
    if CurrentP >= 5:
        CurrentFreeGasViscosity = (mu0 + 1.08e-4 / E * (math.exp(1.44 * r) - math.exp(-1.11 * r ** 1.86))) * 1e-6
    else:
        CurrentFreeGasViscosity = (mu0 + 1.08e-4 / E + (
                0.1023 + 0.23 * r + 0.0585 * r ** 2 - 0.0408 * r ** 3) / E) * 1e-6

    # Плотность жидкости в текущих условиях
    CurrentLiquidDensity = SaturatedOilDensity * (1 - VolumeWater) + CurrentWaterDensity * VolumeWater
    # Поверхностное натяжение нефть - газ
    TensionOilGas = 10 ** -(1.58 + 0.05 * CurrentP) - 72e-6 * (CurrentT - 303)

    # Поверхностное натяжение вода - газ
    TensionWaterGas = 10 ** -(1.19 + 0.01 * CurrentP)

    # Поверхностное натяжение нефть - вода
    TensionOilWater = TensionWaterGas - TensionOilGas

    # Отношение объема свободного газа к объему нефти в условиях трубы
    FreeGasFactor = FreeGasDensity * 1.29 * GasDensity / CurrentFreeGasDensity

    # Плотность смеси Oil+Water+Gas = OWG
    OWG_density = CurrentLiquidDensity + CurrentFreeGasDensity * FreeGasFactor * (1-VolumeWater)

    if not X_kg_sec:
        X_kg_sec = 0
    # Объемный расход нефтеводогазовой смеси в условиях транспорта
    Q_owg = X_kg_sec / OWG_density

    Q_oil = Q_owg * (1 - VolumeWater)
    Q_water = Q_owg * VolumeWater
    Q_liquid = Q_oil + Q_water

    # Объемный расход газовой фазы в условиях транспорта
    Q_gas = Q_oil * FreeGasFactor

    # Объемное расходное газосодержание
    VolumeGas = Q_gas / Q_owg if Q_owg!=0 else 0  # if CurrentP < SaturationPressure_MPa else 0
    Thermal_capacity = VolumeWater * 4183 + (1 - VolumeWater) * 2100
# TODO Separate input and output fluid parameters. It is not necessary to return all, most of them aint used
    return Mishenko(oil_params=calc_params, CurrentP_MPa=CurrentP, CurrentT_K=CurrentT, VolumeWater_fraction=VolumeWater, Q_liq_m3_s=Q_liquid,
                    g=g, SaturationPressure_MPa=Saturation_pressure,
                    DissolvedGasAmount=DissolvedGasAmount, FreeGasDensity_kg_m3=FreeGasDensity,
                    SaturatedOilDensity_kg_m3=SaturatedOilDensity,
                    CurrentWaterDensity_kg_m3=CurrentWaterDensity, CurrentOilViscosity_Pa_s=CurrentOilViscosity,
                    CurrentWaterViscosity_Pa_s=CurrentWaterViscosity,
                    RelativePressure=RelativePressure, RelativeTemperature=RelativeTemperature,
                    CurrentFreeGasDensity_kg_m3=CurrentFreeGasDensity,
                    CurrentFreeGasViscosity_Pa_s=CurrentFreeGasViscosity, CurrentLiquidDensity_kg_m3=CurrentLiquidDensity,
                    TensionOilGas=TensionOilGas,
                    TensionWaterGas=TensionWaterGas, TensionOilWater=TensionOilWater, Q_gas_m3_s=Q_gas, Q_liq_and_gas_m3_s=Q_owg,
                    VolumeGas_fraction=VolumeGas, Thermal_capacity = Thermal_capacity)

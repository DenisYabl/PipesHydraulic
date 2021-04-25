import Tools.HE2_ABC
from Tools.HE2_ABC import oil_params
import GraphEdges.HE2_Pipe as he2_pipe
from GraphEdges.HE2_Plast import HE2_Plast
from GraphEdges.HE2_WellPump import HE2_WellPump
import Fluids.HE2_Fluid as he2_fluid
import GraphNodes.HE2_Vertices

class_aliases = ['OilPipeSegment', 'WaterPipeSegment', 'OilPipe', 'WaterPipe']
class_aliases += ['DummyWater', 'DummyOil', 'OilWater', 'Plast', 'WellPump']

classes = [he2_pipe.HE2_OilPipeSegment, he2_pipe.HE2_WaterPipeSegment, he2_pipe.HE2_OilPipe, he2_pipe.HE2_WaterPipe]
classes += [he2_fluid.HE2_DummyWater, he2_fluid.gimme_dummy_BlackOil, he2_fluid.HE2_BlackOil, HE2_Plast, HE2_WellPump]

get_alias = dict(zip(map(str, classes), class_aliases))
get_class = dict(zip(class_aliases, classes))


def fluid_to_dict(fluid):
    class_alias = get_alias[str(fluid.__class__)]
    if not class_alias in ['DummyWater', 'DummyOil', 'OilWater']:
        return None
    rez = dict(class_alias=class_alias)
    if class_alias in ['DummyOil', 'OilWater']:
        flds = Tools.HE2_ABC.fieldlist
        d = {k: fluid.oil_params.__dict__[k] for k in flds}
        rez.update(d)

    return rez

def dict_to_fluid(obj_dict):
    class_alias = obj_dict.pop('class_alias')

    if class_alias == 'OilWater':
        obj_dict.pop('Q_m3_sec')
        obj_dict.pop('currentP_bar')
        obj_dict.pop('currentT_C')
        obj_dict.pop('CurrentLiquidDensity_kg_m3')
        oil_params = Tools.HE2_ABC.oil_params(**obj_dict)
        rez = he2_fluid.HE2_BlackOil(oil_params)
        return rez

    init_kwargs = dict()
    if class_alias == 'DummyOil':
        init_kwargs = dict(daily_Q = obj_dict['Q_m3_day'], VolumeWater = obj_dict['volumewater_percent'])

    rez = get_class[class_alias](**init_kwargs)
    return rez


def pipesegment_to_dict(pipe_segment):
    class_alias = get_alias[str(pipe_segment.__class__)]
    if not class_alias in ['OilPipeSegment', 'WaterPipeSegment']:
        return None
    rez = dict(class_alias=class_alias)
    rez.update(fluid = fluid_to_dict(pipe_segment.fluid))
    rez.update(inner_diam_m = pipe_segment.inner_diam_m)
    rez.update(roughness_m = pipe_segment.roughness_m)
    rez.update(L_m = pipe_segment.L_m)
    rez.update(uphill_m = pipe_segment.uphill_m)
    return rez

def dict_to_pipesegment(obj_dict):
    fluid = dict_to_fluid(obj_dict.pop('fluid'))
    class_alias = obj_dict.pop('class_alias')
    rez = get_class[class_alias](fluid, **obj_dict)
    return rez


def pipe_to_dict(pipe):
    class_alias = get_alias[str(pipe.__class__)]
    if not class_alias in ['OilPipe', 'WaterPipe']:
        return None
    rez = dict(class_alias=class_alias)
    seg_list = []
    for seg in pipe.segments:
        seg_dict = pipesegment_to_dict(seg)
        seg_list += [seg_dict]
    rez.update(segments=seg_list)
    return rez

def dict_to_pipe(obj_dict):
    seg_list = obj_dict['segments']
    class_alias = obj_dict.pop('class_alias')
    rez = get_class[class_alias]([], [], [], [])
    for seg_dict in seg_list:
        segment = dict_to_pipesegment(seg_dict)
        rez.segments += [segment]
    return rez


def plast_to_dict(plast):
    class_alias = get_alias[str(plast.__class__)]
    if not class_alias in ['Plast']:
        return None
    rez = dict(class_alias=class_alias)
    rez.update(fluid = fluid_to_dict(plast.fluid))
    rez.update(productivity = plast.Productivity)
    return rez

def dict_to_plast(obj_dict):
    class_alias = obj_dict.pop('class_alias')
    productivity = obj_dict.pop('productivity')
    fluid = dict_to_fluid(obj_dict.pop('fluid'))
    rez = get_class[class_alias](productivity, fluid)
    return rez

def wellpump_to_dict(wellpump):
    class_alias = get_alias[str(wellpump.__class__)]
    if not class_alias in ['WellPump']:
        return None
    rez = dict(class_alias=class_alias)


def dict_to_wellpump(obj_dict):
    pass


def test1():
    pseg = he2_pipe.HE2_OilPipeSegment(None, 0.5, 1e-5, 100, 10)
    d = pipesegment_to_dict(pseg)
    pseg2 = dict_to_pipesegment(d)
    print(pseg, pseg2)

def test2():
    op = oil_params.dummy_oil_params(10, 50)
    fluid = he2_fluid.HE2_BlackOil(op)
    pseg = he2_pipe.HE2_OilPipeSegment(fluid, 0.5, 1e-5, 100, 10)
    d = pipesegment_to_dict(pseg)
    pseg2 = dict_to_pipesegment(d)
    print(pseg, pseg2)

def test3():
    op = oil_params.dummy_oil_params(10, 50)
    fluid = he2_fluid.HE2_BlackOil(op)
    pipe = he2_pipe.HE2_OilPipe([100, 200], [-10, 20], [0.3, 0.4], [1e-5, 1e-5], [fluid, fluid])
    d = pipe_to_dict(pipe)
    pipe2 = dict_to_pipe(d)
    print(pipe, pipe2)

def test4():
    op = oil_params.dummy_oil_params(10, 50)
    fluid = he2_fluid.HE2_BlackOil(op)
    plast = HE2_Plast(fluid=fluid, productivity=10)
    d = plast_to_dict(plast)
    plast2 = dict_to_plast(d)
    print(plast, plast2)


if __name__ == '__main__':
    test4()

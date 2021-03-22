import Tools.HE2_ABC as abc
import GraphEdges.HE2_Pipe as he2_pipe
import GraphEdges.HE2_Plast
import GraphEdges.HE2_WellPump
import Fluids.HE2_Fluid as he2_fluid
import GraphNodes.HE2_Vertices

class_aliases = ['OilPipeSegment', 'WaterPipeSegment', 'OilPipe', 'WaterPipe', 'DummyWater', 'DummyOil', 'OilWater']

classes = [he2_pipe.HE2_OilPipeSegment, he2_pipe.HE2_WaterPipeSegment, he2_pipe.HE2_OilPipe, he2_pipe.HE2_WaterPipe]
classes += [he2_fluid.HE2_DummyWater, he2_fluid.HE2_DummyOil, he2_fluid.HE2_OilWater]
get_alias = dict(zip(map(str, classes), class_aliases))
get_class = dict(zip(class_aliases, classes))


def fluid_to_dict(fluid):
    class_alias = get_alias[str(fluid.__class__)]
    if not class_alias in ['DummyWater', 'DummyOil', 'OilWater']:
        return None
    rez = dict(class_alias=class_alias)
    if class_alias in ['DummyOil', 'OilWater']:
        rez.update(fluid.oil_params)
    return rez

def dict_to_fluid(obj_dict):
    class_alias = obj_dict.pop('class_alias')
    rez = get_class[class_alias](**obj_dict)
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

if __name__ == '__main__':
    pseg = he2_pipe.HE2_OilPipeSegment(None, 0.5, 1e-5, 100, 10)
    d = pipesegment_to_dict(pseg)
    pseg2 = dict_to_pipesegment(d)
    print(pseg, pseg2)
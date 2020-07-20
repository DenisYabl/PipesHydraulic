import json
import os
import pandas as pd

def get_global_id():
    if 'ID' not in get_global_id.__dict__:
        get_global_id.ID = 0
    rez = get_global_id.ID
    get_global_id.ID += 1
    return rez


class HE2_Object():
    def __init__(self, obj_id=None, name='', description='', inlets=None, outlets=None, p_link=None):
        if obj_id is None:
            obj_id = get_global_id()
        self.ID = obj_id
        self.Name = name
        self.Description = description
        self.Inlets = dict(Inlet0=None)
        if inlets:
            self.Inlets.update(inlets)
        self.Outlets = dict(Outlet0=None)
        if outlets:
            self.Outlets.update(outlets)

class HE2_Tech_Schema():
    def __init__(self):
        self.objects_by_id = {}
        self.objects_by_full_name = {}
        self.containing_tree = {}
        self.graph = None

    def get_object_by_id(self, id):
        pass

    def get_object_by_full_name(self, id):
        pass

    def get_containing_tree(self, root=None):
        pass

    def get_objects_graph(self):
        pass

    def get_sub_schema(self, id_list=None):
        pass


class HE2_Schema_Persister():
    def __init__(self):
        pass

    def build_tech_shema_from_file(self, filename):
        name, ext = os.path.splitext(filename)
        if ext[0:4] in ('.xls', '.xlsx'):
            lines_colnames = ['line_name', 'node1_name', 'node2_name', 'pipe_num', 'L', 'D', 'Wall', 'Rough']
            df_lines = pd.read_excel(filename, sheet_name=0, header=None, names=lines_colnames, skiprows=range(9))
            print(df_lines)
            pass
        elif ext in ('.txt', '.json'):
            f = open(filename, 'r', encoding='UTF-8')
            ts_json = json.load(f)
            rez = self.build_tech_shema_from_json(ts_json)
            return rez

    def build_tech_shema_from_json(self, ts_json):
        pass

    def build_tech_shema_from_dataframes(self, nodes_df, pipes_df):
        pass

    def dump_tech_schema_to_json(self, tech_schema):
        pass

    def dump_tech_schema_to_dataframes(self, tech_schema):
        pass


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    # def_fn = r'..\data\tech_schema.json'
    def_fn = r'..\data\waterpipes.xlsx'
    parser.add_argument("-f", "--filename", type=str, required=False, help="Techshema json file", default=def_fn)
    args = parser.parse_args()
    filename = args.filename
    persist = HE2_Schema_Persister()
    tech_schema = persist.build_tech_shema_from_file(filename)
    print(tech_schema)

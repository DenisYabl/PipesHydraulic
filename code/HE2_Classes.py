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

    def build_tech_shema_from_json(self, ts_json):
        pass

    def build_tech_shema_from_dataframes(self, nodes_df, pipes_df):
        pass

    def dump_tech_schema_to_json(self, tech_schema):
        pass

    def dump_tech_schema_to_dataframes(self, tech_schema):
        pass


if __name__ == '__main__':
    hello_world = HE2_Object()
    print(hello_world.ID)


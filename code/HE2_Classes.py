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


if __name__ == '__main__':
    hello_world = HE2_Object()
    print(hello_world.ID)


import pandas as pd

def load_data_from_dataset(dataset, schema):
    """
    Convert dataset to pandas.DataFrame
    Args:
        dataset: list of lists of dataset rows
        schema: schema of the dataset

    Returns:
    """
    idxs = [i for i in range(1, len(schema.split('`'))) if i % 2 != 0]
    column_names = [schema.split("`")[i] for i in idxs]
    df = pd.DataFrame(dataset, columns=column_names)
    return df

def load_dataset_and_schema(dataFrame):
    dataset = []
    columns = dataFrame.columns
    for _, rows in dataFrame.iterrows():
        row = []
        for col in columns:
           row.append(rows[col])
        dataset.append(row)
    schema = ""
    for col in columns:
        schema += f"`{col}` {str(dataFrame[col].dtype).upper()}, "
    return dataset, schema

def load_well_params(data_path):
    df = pd.read_csv(data_path)
    return df
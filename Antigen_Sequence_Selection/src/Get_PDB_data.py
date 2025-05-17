from rcsbapi.data import DataSchema, DataQuery

from rcsbapi.data import DataQuery as Query

query = Query(
    input_type="entries",
    input_ids=["4HHB"],
    return_data_list=["exptl.method"]
)

result_dict = query.exec()
print(result_dict)
# print(query.get_response()) would be equivalent
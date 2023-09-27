import numpy as np

mols=np.array([[[range(0),
         range(1),
         range(2)],
        [range(3),
         range(4),
         range(5)]],

       [[range(6),
         range(7),
         range(8)],
        [range(9),
         range(10),
         range(11)]]], dtype=object)

for block in mols:
    for row in block:
        for mol in row:
            print(mol)

def print_ndarray_elements(arr, prefix=''):
    if len(arr.shape) == 0:
        print(prefix + str(arr))
    else:
        for i in range(arr.shape[0]):
            new_prefix = prefix + "[" + str(i) + "]"
            print_ndarray_elements(arr[i], new_prefix)

print_ndarray_elements(mols)

# for item in mols:
#     current_level_item = item
#     if isinstance(current_level_item, np.ndarray):
#         while isinstance(current_level_item, np.ndarray):
#                     for current_level_index, current_level_item in enumerate(current_level_item):
#             print(f"{level=} {current_level_index=}, {current_level_item=}, is ndarray={isinstance(item, np.ndarray)}")


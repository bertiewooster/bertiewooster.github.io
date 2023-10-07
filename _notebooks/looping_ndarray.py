import numpy as np


# Define the simplified MolecularIsotope class with only the 'name' property
class MolecularIsotope:
    def __init__(self, name):
        self.name = name

# Define the number of elements you want in the array
num_elements = 5  # You can change this to any desired number

# Create a 1-dimensional ndarray with the specified number of elements and element type of object
arr = np.empty(num_elements, dtype=object)

# Create separate instances of the simplified MolecularIsotope class for each element
for i in range(num_elements):
    arr[i] = MolecularIsotope("Hydrogen")

# Print the elements of the ndarray and their memory addresses
for i, isotope in enumerate(arr):
    print(f"Element {i}: {isotope.name}")
    print(f"Memory Address: {id(isotope)}")

interleaved = [['quadrilateral', '', '', ''], ['|', '|', '', '|'], ['trapezoid', 'parallelogram', '', 'kite'], ['|', '|', '|', ''], ['isosceles trapezoid', 'rhombus', 'rectangle', ''], ['', '|', '|', ''], ['', 'square', 'square', '']]

f_string = ''
for col_index, col in enumerate(interleaved[0]):
    f_string += "{row[" + str(col_index) + "]:<20}"

new_line = '\n'
tree_plot = ""

for row in interleaved:
    # print(f_string.format(row=row))
    tree_plot += f_string.format(row=row) + new_line

print(tree_plot)
    # return tree_plot

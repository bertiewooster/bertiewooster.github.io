```python
import matplotlib.pyplot as plt
```


```python
a = 7
```


```python
b = a + 3
b
```




    10




```python
# Create the plot
fig, ax = plt.subplots()
ax.scatter(a, b)
ax.set_xlabel("a")
ax.set_ylabel("b")
ax.set_title("Plot of b against a")
```




    Text(0.5, 1.0, 'Plot of b against a')




    
![png](/images/2026-01-05-ChEBML-database_files/2026-01-05-ChEBML-database_3_1.png)
    


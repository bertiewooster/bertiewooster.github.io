# Chemistry machine learning for drug discovery with DeepChem

This example uses machine learning to predict the lipophilicty of compounds.

[Lipophilicty](https://en.wikipedia.org/wiki/Lipophilicity) measures how well a compound dissolves in non-polar media such as fats and lipids. So it's important for drugs that are delivered orally (for example, via a pill) because the active ingredient [needs to be absorbed into the lipids](https://emerypharma.com/blog/drug-lipophilicity-and-absorption-a-continuous-challenge-toward-the-goal-of-drug-discovery/) of biological membranes.

_[Open this notebook in Google Colab](https://colab.research.google.com/drive/1P7txNo9M6cln-iukWRrLkqFFrZvx7jI3?usp=sharing)_


```python
%%capture
!pip install --pre deepchem[tensorflow]
```


```python
#For mol_frame
%%capture
!pip install git+https://github.com/apahl/mol_frame

from mol_frame import mol_frame as mf
import os, holoviews as hv
os.environ['HV_DOC_HTML'] = 'true'
```


```python
import deepchem as dc
import seaborn
import pandas as pd
import matplotlib.pyplot as plt
```

[DeepChem](https://deepchem.io/) is a free and open-source Python package for deep learning for chemistry and other sciences. DeepChem has a [lipophilicty data set](https://deepchem.readthedocs.io/en/latest/api_reference/moleculenet.html#lipo-datasets) contains measured [logD](https://www.cambridgemedchemconsulting.com/resources/physiochem/logD.html) values for 4200 compounds.

As usual for machine learning (ML), we split the data set into training and test data. We train the ML model on the train data, then apply it to the test data and check how well the model predicts the lipophilicty of compounds that the model hasn't processed before.

For this data set, we [split by scaffold](https://deepchem.readthedocs.io/en/latest/api_reference/splitters.html#scaffoldsplitter) the 4200 compounds based on the [Bemis-Murcko scaffold representation](https://pubs.acs.org/doi/10.1021/jm9602928). Such splitting groups molecules based on their scaffolds (core structure) to [prevent train and test from having very similar molecules](httphttps://github.com/deepchem/deepchem/blob/master/examples/tutorials/Working_With_Splitters.ipynb), which could lead to the model appearing to perform well on the test set, but then performing poorly on less-similar molecules in production.


```python
tasks, datasets, transformers = dc.molnet.load_lipo(featurizer='GraphConv', splitter='Scaffold')
train_dataset, valid_dataset, test_dataset = datasets
```

The number of compounds in the train, validate, and test sets is:


```python
splits = (train_dataset.X.shape[0], valid_dataset.X.shape[0], test_dataset.X.shape[0])
splits
```




    (3360, 420, 420)



Which represents an 80:10:10 train:validate:test split:


```python
[split / sum(splits) for split in splits]
```




    [0.8, 0.1, 0.1]



Next, we build a model using DeepChem's [graph convolutional network](https://deepchem.readthedocs.io/en/latest/api_reference/models.html). We use the [dropout technique](https://ai-pool.com/a/s/dropout-in-deep-learning) to avoid overfitting.


```python
model = dc.models.GraphConvModel(n_tasks=1, mode='regression', dropout=0.2)
```

Then we train the model on the train dataset.


```python
# Note that training with 200 epochs takes around 8 minutes on Google Colab.
# If you are debugging or just running a proof of concept, you may want to set nb_epoch=10 to speed execution.
%%capture
model.fit(train_dataset, nb_epoch=200);
```

To check how well the model fits the train and test data, we examine the [Pearson correlation coefficient](https://www.scribbr.com/statistics/pearson-correlation-coefficient/) score. Its value can range from -1 to 1, where positive values indicate postive correlation, zero indicates no correlation, and negative values indicate negative correlation. The magnitude of the value indicates the strength of the correlation: less than 0.3 is weak, 0.3-0.5 is moderate, and 0.5-1 is strong.


```python
metric = dc.metrics.Metric(dc.metrics.pearson_r2_score)
print("Training set score:", model.evaluate(train_dataset, [metric], transformers))
print("Test set score:", model.evaluate(test_dataset, [metric], transformers))
```

    Training set score: {'pearson_r2_score': 0.8866760139249386}
    Test set score: {'pearson_r2_score': 0.5358381141161308}


The [Pearson correlation coefficient](https://www.scribbr.com/statistics/pearson-correlation-coefficient/) score is worse for the test data than for the train data. This is expected because the test data is new to the model. Nevertheless, the magnitude being greater than 0.5 indicates a strong, positive correlation betwen predicted and measued lipophilicity.

One contributing factor might be that our train set may not have molecules similar enough to those in the test set. Recall that we [split by scaffold](https://deepchem.readthedocs.io/en/latest/api_reference/splitters.html#scaffoldsplitter), so it's possible that such splitting led to compounds in the test set that have scaffolds significantly different from those in the train set.

Adding compounds to the dataset so that there is less "scaffold distance" (difference in scaffold structure) between groups of compounds might help. To take a simple hypothetical example, if we had a dataset with compounds containing fused rings, with many two-ring compounds and few four-ring compounds, scaffold splitting might put all the two-ring compounds in the train set and all the four-ring compounds in the test set. We expect that would lead to poor predictions on the test set. We would want to augment the dataset by adding compounds containing three fused rings.

To learn more about how well the model works, let's apply it to the test data and compare the predicted and experimental results.


```python
lipos = model.predict_on_batch(test_dataset.X)
```

Then we make put the measured and model-predicted results into a [pandas dataframe](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) for easy processing.


```python
lipo_list = []
expt_lipo_list = []
for molecule, lipo, test_lipo in zip(test_dataset.ids, lipos, test_dataset.y):
    lipo_list += [lipo[0]]
    expt_lipo_list += [test_lipo[0]]
    
compound_id = list(range(len(test_dataset.ids)))
df = pd.DataFrame(list(zip(expt_lipo_list, lipo_list, test_dataset.ids, compound_id)), columns = ["measured", "predicted", "Smiles", "Compound_Id"])
df
```





  <div id="df-c1bdccd2-9722-4878-9a57-d79026309bf5">
    <div class="colab-df-container">
      <div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>measured</th>
      <th>predicted</th>
      <th>Smiles</th>
      <th>Compound_Id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>-1.810832</td>
      <td>-0.892024</td>
      <td>O[C@@H](CNCCCOCCNCCc1cccc(Cl)c1)c2ccc(O)c3NC(=O)Sc23</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.319651</td>
      <td>-0.728260</td>
      <td>NC1=NN(CC1)c2cccc(c2)C(F)(F)F</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>-0.192325</td>
      <td>0.023027</td>
      <td>COc1cc(ccc1Cn2ccc3ccc(cc23)C(=O)NCC4CCCC4)C(=O)NS(=O)(=O)c5ccccc5</td>
      <td>2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.938978</td>
      <td>0.337164</td>
      <td>NC(=O)NC(=O)C(Nc1ccc2CCCc2c1)c3ccccc3</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.856401</td>
      <td>0.464554</td>
      <td>COc1ccc(cn1)c2ccc3nc(N)sc3c2</td>
      <td>4</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>415</th>
      <td>0.815112</td>
      <td>1.280881</td>
      <td>COc1cc(ccc1N2CC[C@@H](O)C2)N3N=Nc4cc(sc4C3=O)c5ccc(Cl)cc5</td>
      <td>415</td>
    </tr>
    <tr>
      <th>416</th>
      <td>1.327089</td>
      <td>0.092211</td>
      <td>CCC(COC(=O)c1cc(OC)c(OC)c(OC)c1)(N(C)C)c2ccccc2</td>
      <td>416</td>
    </tr>
    <tr>
      <th>417</th>
      <td>-0.175810</td>
      <td>0.370334</td>
      <td>CCCSc1ncccc1C(=O)N2CCCC2c3ccncc3</td>
      <td>417</td>
    </tr>
    <tr>
      <th>418</th>
      <td>0.071921</td>
      <td>-0.493196</td>
      <td>Oc1ncnc2scc(c3ccsc3)c12</td>
      <td>418</td>
    </tr>
    <tr>
      <th>419</th>
      <td>0.806855</td>
      <td>0.264380</td>
      <td>OC1(CN2CCC1CC2)C#Cc3ccc(cc3)c4ccccc4</td>
      <td>419</td>
    </tr>
  </tbody>
</table>
<p>420 rows × 4 columns</p>
</div>
      <button class="colab-df-convert" onclick="convertToInteractive('df-c1bdccd2-9722-4878-9a57-d79026309bf5')"
              title="Convert this dataframe to an interactive table."
              style="display:none;">

  <svg xmlns="http://www.w3.org/2000/svg" height="24px"viewBox="0 0 24 24"
       width="24px">
    <path d="M0 0h24v24H0V0z" fill="none"/>
    <path d="M18.56 5.44l.94 2.06.94-2.06 2.06-.94-2.06-.94-.94-2.06-.94 2.06-2.06.94zm-11 1L8.5 8.5l.94-2.06 2.06-.94-2.06-.94L8.5 2.5l-.94 2.06-2.06.94zm10 10l.94 2.06.94-2.06 2.06-.94-2.06-.94-.94-2.06-.94 2.06-2.06.94z"/><path d="M17.41 7.96l-1.37-1.37c-.4-.4-.92-.59-1.43-.59-.52 0-1.04.2-1.43.59L10.3 9.45l-7.72 7.72c-.78.78-.78 2.05 0 2.83L4 21.41c.39.39.9.59 1.41.59.51 0 1.02-.2 1.41-.59l7.78-7.78 2.81-2.81c.8-.78.8-2.07 0-2.86zM5.41 20L4 18.59l7.72-7.72 1.47 1.35L5.41 20z"/>
  </svg>
      </button>

  <style>
    .colab-df-container {
      display:flex;
      flex-wrap:wrap;
      gap: 12px;
    }

    .colab-df-convert {
      background-color: #E8F0FE;
      border: none;
      border-radius: 50%;
      cursor: pointer;
      display: none;
      fill: #1967D2;
      height: 32px;
      padding: 0 0 0 0;
      width: 32px;
    }

    .colab-df-convert:hover {
      background-color: #E2EBFA;
      box-shadow: 0px 1px 2px rgba(60, 64, 67, 0.3), 0px 1px 3px 1px rgba(60, 64, 67, 0.15);
      fill: #174EA6;
    }

    [theme=dark] .colab-df-convert {
      background-color: #3B4455;
      fill: #D2E3FC;
    }

    [theme=dark] .colab-df-convert:hover {
      background-color: #434B5C;
      box-shadow: 0px 1px 3px 1px rgba(0, 0, 0, 0.15);
      filter: drop-shadow(0px 1px 2px rgba(0, 0, 0, 0.3));
      fill: #FFFFFF;
    }
  </style>

      <script>
        const buttonEl =
          document.querySelector('#df-c1bdccd2-9722-4878-9a57-d79026309bf5 button.colab-df-convert');
        buttonEl.style.display =
          google.colab.kernel.accessAllowed ? 'block' : 'none';

        async function convertToInteractive(key) {
          const element = document.querySelector('#df-c1bdccd2-9722-4878-9a57-d79026309bf5');
          const dataTable =
            await google.colab.kernel.invokeFunction('convertToInteractive',
                                                     [key], {});
          if (!dataTable) return;

          const docLinkHtml = 'Like what you see? Visit the ' +
            '<a target="_blank" href=https://colab.research.google.com/notebooks/data_table.ipynb>data table notebook</a>'
            + ' to learn more about interactive tables.';
          element.innerHTML = '';
          dataTable['output_type'] = 'display_data';
          await google.colab.output.renderOutput(dataTable, element);
          const docLink = document.createElement('div');
          docLink.innerHTML = docLinkHtml;
          element.appendChild(docLink);
        }
      </script>
    </div>
  </div>




Now we can use a scatter plot to compare the predicted against measured values. We use the [seaborn statistical data visualization package](https://seaborn.pydata.org/) to plot the data. We show the line where the predicted and measured lipophilicty values are equal, in other words the line that all points would lie on if the model made perfect predictions.


```python
seaborn.scatterplot(data=df, x = "measured", y = "predicted").set(title='Lipophilicty for test data:\noctanol/water distribution coefficient\n(logD at pH 7.4)\n');
equal_line_x = [-3, 2]
equal_line_y = equal_line_x
plt.plot(equal_line_x, equal_line_y, color='k', linewidth=0.5);
```


    
![Predicted against measured lipophilicty for test data](/images/2022-12-13-Chemistry-machine-learning-for-drug-discovery-with-DeepChem_files/2022-12-13-Chemistry-machine-learning-for-drug-discovery-with-DeepChem_24_0.png)
    


We can also use [mol_frame](https://github.com/apahl/mol_frame) to make an interactive plot: If you mouse over a marker on the graph, the molecular structure will be displayed! (One interesting thing to note is that molecules which are close together on the plot may have very different molecular structures.) Also, you can pan and zoom the graph, and save it.

*Unfortunately, the interactive plot is not working in the blog version of this notebook. Please visit the [Google Colab notebook](https://colab.research.google.com/drive/1P7txNo9M6cln-iukWRrLkqFFrZvx7jI3?usp=sharing) to access the interactive plot.*


```python
cpds = mf.MolFrame(df)
cpds = cpds.add_b64()
hv.extension('bokeh')
cpds.scatter("measured", "predicted")
```

    * using Smiles
    * add b64:               (  420 |    5)

    * using Mol_b64
    * add img:               (  420 |    6)


Let's plot the train data set to visually compare those two sets. We made the markers (points) smaller because there are so many and we don't want them to overlap with each other so much.


```python
lipos_train = model.predict_on_batch(train_dataset.X)
lipo_list_train = []
expt_lipo_list_train = []
for molecule, lipo, test_lipo in zip(train_dataset.ids, lipos_train, train_dataset.y):
    lipo_list_train += [lipo[0]]
    expt_lipo_list_train += [test_lipo[0]]
df_train = pd.DataFrame(list(zip(expt_lipo_list_train, lipo_list_train)), columns = ["measured", "predicted"])
seaborn.scatterplot(data=df_train, x = "measured", y = "predicted", s=5, palette=["orange"]).set(title='Lipophilicty for train data:\noctanol/water distribution coefficient\n(logD at pH 7.4)\n');
plt.plot(equal_line_x, equal_line_y, color='k', linewidth=0.5);
```


    
![Predicted against measured lipophilicty for train data](/images/2022-12-13-Chemistry-machine-learning-for-drug-discovery-with-DeepChem_files/2022-12-13-Chemistry-machine-learning-for-drug-discovery-with-DeepChem_28_0.png)
    


*Unfortunately, mol_frame takes too long to render thousands of molecules so that Google Colab disconnects before it completes. So we're not showing the interactive plot here. You probably could if you ran this notebook locally.*

To overlay the two plots, we'll first [concatenate the data](https://stackoverflow.com/questions/51732867/seaborn-plot-two-data-sets-on-the-same-scatter-plot#51733133) in the two sets (test and train), adding a `dataset` column to keep track of which set each point came from.


```python
# Remove the columns we added to the test dataframe for mol_frame
try:
  df.drop(columns=["Smiles", "Compound_Id"], inplace=True)
except:
  # Suppress any error if those columns have already been dropped
  pass

concatenated = pd.concat([df.assign(dataset='test'), df_train.assign(dataset='train')])
concatenated
```





  <div id="df-236d7f19-6b5b-47ac-b707-7d324b8e28ef">
    <div class="colab-df-container">
      <div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>measured</th>
      <th>predicted</th>
      <th>dataset</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>-1.810832</td>
      <td>-0.892024</td>
      <td>test</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.319651</td>
      <td>-0.728260</td>
      <td>test</td>
    </tr>
    <tr>
      <th>2</th>
      <td>-0.192325</td>
      <td>0.023027</td>
      <td>test</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.938978</td>
      <td>0.337164</td>
      <td>test</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.856401</td>
      <td>0.464554</td>
      <td>test</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>3355</th>
      <td>0.526093</td>
      <td>-0.117241</td>
      <td>train</td>
    </tr>
    <tr>
      <th>3356</th>
      <td>-0.027172</td>
      <td>-0.238245</td>
      <td>train</td>
    </tr>
    <tr>
      <th>3357</th>
      <td>-1.629163</td>
      <td>-1.749196</td>
      <td>train</td>
    </tr>
    <tr>
      <th>3358</th>
      <td>0.633443</td>
      <td>0.961528</td>
      <td>train</td>
    </tr>
    <tr>
      <th>3359</th>
      <td>-0.960290</td>
      <td>-1.251761</td>
      <td>train</td>
    </tr>
  </tbody>
</table>
<p>3780 rows × 3 columns</p>
</div>
      <button class="colab-df-convert" onclick="convertToInteractive('df-236d7f19-6b5b-47ac-b707-7d324b8e28ef')"
              title="Convert this dataframe to an interactive table."
              style="display:none;">

  <svg xmlns="http://www.w3.org/2000/svg" height="24px"viewBox="0 0 24 24"
       width="24px">
    <path d="M0 0h24v24H0V0z" fill="none"/>
    <path d="M18.56 5.44l.94 2.06.94-2.06 2.06-.94-2.06-.94-.94-2.06-.94 2.06-2.06.94zm-11 1L8.5 8.5l.94-2.06 2.06-.94-2.06-.94L8.5 2.5l-.94 2.06-2.06.94zm10 10l.94 2.06.94-2.06 2.06-.94-2.06-.94-.94-2.06-.94 2.06-2.06.94z"/><path d="M17.41 7.96l-1.37-1.37c-.4-.4-.92-.59-1.43-.59-.52 0-1.04.2-1.43.59L10.3 9.45l-7.72 7.72c-.78.78-.78 2.05 0 2.83L4 21.41c.39.39.9.59 1.41.59.51 0 1.02-.2 1.41-.59l7.78-7.78 2.81-2.81c.8-.78.8-2.07 0-2.86zM5.41 20L4 18.59l7.72-7.72 1.47 1.35L5.41 20z"/>
  </svg>
      </button>

  <style>
    .colab-df-container {
      display:flex;
      flex-wrap:wrap;
      gap: 12px;
    }

    .colab-df-convert {
      background-color: #E8F0FE;
      border: none;
      border-radius: 50%;
      cursor: pointer;
      display: none;
      fill: #1967D2;
      height: 32px;
      padding: 0 0 0 0;
      width: 32px;
    }

    .colab-df-convert:hover {
      background-color: #E2EBFA;
      box-shadow: 0px 1px 2px rgba(60, 64, 67, 0.3), 0px 1px 3px 1px rgba(60, 64, 67, 0.15);
      fill: #174EA6;
    }

    [theme=dark] .colab-df-convert {
      background-color: #3B4455;
      fill: #D2E3FC;
    }

    [theme=dark] .colab-df-convert:hover {
      background-color: #434B5C;
      box-shadow: 0px 1px 3px 1px rgba(0, 0, 0, 0.15);
      filter: drop-shadow(0px 1px 2px rgba(0, 0, 0, 0.3));
      fill: #FFFFFF;
    }
  </style>

      <script>
        const buttonEl =
          document.querySelector('#df-236d7f19-6b5b-47ac-b707-7d324b8e28ef button.colab-df-convert');
        buttonEl.style.display =
          google.colab.kernel.accessAllowed ? 'block' : 'none';

        async function convertToInteractive(key) {
          const element = document.querySelector('#df-236d7f19-6b5b-47ac-b707-7d324b8e28ef');
          const dataTable =
            await google.colab.kernel.invokeFunction('convertToInteractive',
                                                     [key], {});
          if (!dataTable) return;

          const docLinkHtml = 'Like what you see? Visit the ' +
            '<a target="_blank" href=https://colab.research.google.com/notebooks/data_table.ipynb>data table notebook</a>'
            + ' to learn more about interactive tables.';
          element.innerHTML = '';
          dataTable['output_type'] = 'display_data';
          await google.colab.output.renderOutput(dataTable, element);
          const docLink = document.createElement('div');
          docLink.innerHTML = docLinkHtml;
          element.appendChild(docLink);
        }
      </script>
    </div>
  </div>




Now we can plot the two datasets on one graph.


```python
seaborn.scatterplot(x='measured', 
                    y='predicted', 
                    data=concatenated,
                    hue='dataset',
                    s=5);
plt.plot(equal_line_x, equal_line_y, color='k', linewidth=0.5);
```


    
![Predicted against measured lipophilicty for test and train data](/images/2022-12-13-Chemistry-machine-learning-for-drug-discovery-with-DeepChem_files/2022-12-13-Chemistry-machine-learning-for-drug-discovery-with-DeepChem_33_0.png)
    


Some of the test (blue) data points are predicted (on the vertical axis) outliers, reflecting that the model performs poorly for them. We might consider featurizing our data further to let the model predict lipophilicity based on more properties of the compounds. Above, we used the [`GraphConv`](https://deepchem.readthedocs.io/en/latest/api_reference/featurizers.html#graph-convolution-featurizers) featurizer, which represents only the atoms in a molecule. We might try the [`WeaveFeaturizer`](https://deepchem.readthedocs.io/en/latest/api_reference/featurizers.html#weavefeaturizer) which also represents the bonds, though it requires more resources because it stores the relationship between each pair of atoms in a molecule.

This blog post was based on the DeepChem tutorial [The Basic Tools of the Deep Life Sciences](https://github.com/deepchem/deepchem/blob/master/examples/tutorials/The_Basic_Tools_of_the_Deep_Life_Sciences.ipynb). Thanks to [mol_frame](https://github.com/apahl/mol_frame) author [Axel Pahl](https://github.com/apahl) for making that package usable from Google Colab.

This post was updated on December 28, 2022 to explain the [Pearson correlation coefficient](https://www.scribbr.com/statistics/pearson-correlation-coefficient/) and increase the number of training epochs from 100 to 200.

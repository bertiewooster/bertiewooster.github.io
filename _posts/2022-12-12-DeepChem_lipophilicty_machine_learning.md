```python
!pip install --pre deepchem[tensorflow]
```

    Looking in indexes: https://pypi.org/simple, https://us-python.pkg.dev/colab-wheels/public/simple/
    Collecting deepchem[tensorflow]
      Downloading deepchem-2.7.1-py3-none-any.whl (693 kB)
    [K     |████████████████████████████████| 693 kB 6.0 MB/s 
    [?25hCollecting rdkit
      Downloading rdkit-2022.9.3-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (29.3 MB)
    [K     |████████████████████████████████| 29.3 MB 6.8 MB/s 
    [?25hRequirement already satisfied: joblib in /usr/local/lib/python3.8/dist-packages (from deepchem[tensorflow]) (1.2.0)
    Requirement already satisfied: pandas in /usr/local/lib/python3.8/dist-packages (from deepchem[tensorflow]) (1.3.5)
    Requirement already satisfied: numpy>=1.21 in /usr/local/lib/python3.8/dist-packages (from deepchem[tensorflow]) (1.21.6)
    Requirement already satisfied: scipy<1.9 in /usr/local/lib/python3.8/dist-packages (from deepchem[tensorflow]) (1.7.3)
    Requirement already satisfied: scikit-learn in /usr/local/lib/python3.8/dist-packages (from deepchem[tensorflow]) (1.0.2)
    Requirement already satisfied: tensorflow in /usr/local/lib/python3.8/dist-packages (from deepchem[tensorflow]) (2.9.2)
    Requirement already satisfied: tensorflow-probability in /usr/local/lib/python3.8/dist-packages (from deepchem[tensorflow]) (0.17.0)
    Collecting tensorflow-addons
      Downloading tensorflow_addons-0.19.0-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl (1.1 MB)
    [K     |████████████████████████████████| 1.1 MB 27.5 MB/s 
    [?25hRequirement already satisfied: python-dateutil>=2.7.3 in /usr/local/lib/python3.8/dist-packages (from pandas->deepchem[tensorflow]) (2.8.2)
    Requirement already satisfied: pytz>=2017.3 in /usr/local/lib/python3.8/dist-packages (from pandas->deepchem[tensorflow]) (2022.6)
    Requirement already satisfied: six>=1.5 in /usr/local/lib/python3.8/dist-packages (from python-dateutil>=2.7.3->pandas->deepchem[tensorflow]) (1.15.0)
    Requirement already satisfied: Pillow in /usr/local/lib/python3.8/dist-packages (from rdkit->deepchem[tensorflow]) (7.1.2)
    Requirement already satisfied: threadpoolctl>=2.0.0 in /usr/local/lib/python3.8/dist-packages (from scikit-learn->deepchem[tensorflow]) (3.1.0)
    Requirement already satisfied: keras<2.10.0,>=2.9.0rc0 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (2.9.0)
    Requirement already satisfied: libclang>=13.0.0 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (14.0.6)
    Requirement already satisfied: wrapt>=1.11.0 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (1.14.1)
    Requirement already satisfied: setuptools in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (57.4.0)
    Requirement already satisfied: termcolor>=1.1.0 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (2.1.1)
    Requirement already satisfied: keras-preprocessing>=1.1.1 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (1.1.2)
    Requirement already satisfied: opt-einsum>=2.3.2 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (3.3.0)
    Requirement already satisfied: google-pasta>=0.1.1 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (0.2.0)
    Requirement already satisfied: flatbuffers<2,>=1.12 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (1.12)
    Requirement already satisfied: grpcio<2.0,>=1.24.3 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (1.51.1)
    Requirement already satisfied: protobuf<3.20,>=3.9.2 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (3.19.6)
    Requirement already satisfied: typing-extensions>=3.6.6 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (4.4.0)
    Requirement already satisfied: tensorflow-estimator<2.10.0,>=2.9.0rc0 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (2.9.0)
    Requirement already satisfied: tensorflow-io-gcs-filesystem>=0.23.1 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (0.28.0)
    Requirement already satisfied: packaging in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (21.3)
    Requirement already satisfied: h5py>=2.9.0 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (3.1.0)
    Requirement already satisfied: absl-py>=1.0.0 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (1.3.0)
    Requirement already satisfied: gast<=0.4.0,>=0.2.1 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (0.4.0)
    Requirement already satisfied: astunparse>=1.6.0 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (1.6.3)
    Requirement already satisfied: tensorboard<2.10,>=2.9 in /usr/local/lib/python3.8/dist-packages (from tensorflow->deepchem[tensorflow]) (2.9.1)
    Requirement already satisfied: wheel<1.0,>=0.23.0 in /usr/local/lib/python3.8/dist-packages (from astunparse>=1.6.0->tensorflow->deepchem[tensorflow]) (0.38.4)
    Requirement already satisfied: werkzeug>=1.0.1 in /usr/local/lib/python3.8/dist-packages (from tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (1.0.1)
    Requirement already satisfied: google-auth<3,>=1.6.3 in /usr/local/lib/python3.8/dist-packages (from tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (2.15.0)
    Requirement already satisfied: markdown>=2.6.8 in /usr/local/lib/python3.8/dist-packages (from tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (3.4.1)
    Requirement already satisfied: tensorboard-data-server<0.7.0,>=0.6.0 in /usr/local/lib/python3.8/dist-packages (from tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (0.6.1)
    Requirement already satisfied: google-auth-oauthlib<0.5,>=0.4.1 in /usr/local/lib/python3.8/dist-packages (from tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (0.4.6)
    Requirement already satisfied: tensorboard-plugin-wit>=1.6.0 in /usr/local/lib/python3.8/dist-packages (from tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (1.8.1)
    Requirement already satisfied: requests<3,>=2.21.0 in /usr/local/lib/python3.8/dist-packages (from tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (2.23.0)
    Requirement already satisfied: pyasn1-modules>=0.2.1 in /usr/local/lib/python3.8/dist-packages (from google-auth<3,>=1.6.3->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (0.2.8)
    Requirement already satisfied: rsa<5,>=3.1.4 in /usr/local/lib/python3.8/dist-packages (from google-auth<3,>=1.6.3->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (4.9)
    Requirement already satisfied: cachetools<6.0,>=2.0.0 in /usr/local/lib/python3.8/dist-packages (from google-auth<3,>=1.6.3->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (5.2.0)
    Requirement already satisfied: requests-oauthlib>=0.7.0 in /usr/local/lib/python3.8/dist-packages (from google-auth-oauthlib<0.5,>=0.4.1->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (1.3.1)
    Requirement already satisfied: importlib-metadata>=4.4 in /usr/local/lib/python3.8/dist-packages (from markdown>=2.6.8->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (4.13.0)
    Requirement already satisfied: zipp>=0.5 in /usr/local/lib/python3.8/dist-packages (from importlib-metadata>=4.4->markdown>=2.6.8->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (3.11.0)
    Requirement already satisfied: pyasn1<0.5.0,>=0.4.6 in /usr/local/lib/python3.8/dist-packages (from pyasn1-modules>=0.2.1->google-auth<3,>=1.6.3->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (0.4.8)
    Requirement already satisfied: chardet<4,>=3.0.2 in /usr/local/lib/python3.8/dist-packages (from requests<3,>=2.21.0->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (3.0.4)
    Requirement already satisfied: urllib3!=1.25.0,!=1.25.1,<1.26,>=1.21.1 in /usr/local/lib/python3.8/dist-packages (from requests<3,>=2.21.0->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (1.24.3)
    Requirement already satisfied: certifi>=2017.4.17 in /usr/local/lib/python3.8/dist-packages (from requests<3,>=2.21.0->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (2022.9.24)
    Requirement already satisfied: idna<3,>=2.5 in /usr/local/lib/python3.8/dist-packages (from requests<3,>=2.21.0->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (2.10)
    Requirement already satisfied: oauthlib>=3.0.0 in /usr/local/lib/python3.8/dist-packages (from requests-oauthlib>=0.7.0->google-auth-oauthlib<0.5,>=0.4.1->tensorboard<2.10,>=2.9->tensorflow->deepchem[tensorflow]) (3.2.2)
    Requirement already satisfied: pyparsing!=3.0.5,>=2.0.2 in /usr/local/lib/python3.8/dist-packages (from packaging->tensorflow->deepchem[tensorflow]) (3.0.9)
    Requirement already satisfied: typeguard>=2.7 in /usr/local/lib/python3.8/dist-packages (from tensorflow-addons->deepchem[tensorflow]) (2.7.1)
    Requirement already satisfied: cloudpickle>=1.3 in /usr/local/lib/python3.8/dist-packages (from tensorflow-probability->deepchem[tensorflow]) (1.5.0)
    Requirement already satisfied: decorator in /usr/local/lib/python3.8/dist-packages (from tensorflow-probability->deepchem[tensorflow]) (4.4.2)
    Requirement already satisfied: dm-tree in /usr/local/lib/python3.8/dist-packages (from tensorflow-probability->deepchem[tensorflow]) (0.1.7)
    Installing collected packages: rdkit, tensorflow-addons, deepchem
    Successfully installed deepchem-2.7.1 rdkit-2022.9.3 tensorflow-addons-0.19.0



```python
import deepchem as dc
dc.__version__
```

    WARNING:deepchem.models.torch_models:Skipped loading modules with pytorch-geometric dependency, missing a dependency. No module named 'torch_geometric'
    WARNING:deepchem.models:Skipped loading modules with pytorch-geometric dependency, missing a dependency. cannot import name 'DMPNN' from 'deepchem.models.torch_models' (/usr/local/lib/python3.8/dist-packages/deepchem/models/torch_models/__init__.py)
    WARNING:deepchem.models:Skipped loading modules with pytorch-lightning dependency, missing a dependency. No module named 'pytorch_lightning'
    WARNING:deepchem.models:Skipped loading some Jax models, missing a dependency. No module named 'haiku'





    '2.7.1'



DeepChem's [lipophilicty data set](https://deepchem.readthedocs.io/en/latest/api_reference/moleculenet.html#lipo-datasets) contains measured [logD](https://www.cambridgemedchemconsulting.com/resources/physiochem/logD.html) values for 4200 compounds. [Lipophilicty](https://en.wikipedia.org/wiki/Lipophilicity) measures how well a compound dissolves in non-polar media such as fats and lipids. So it's important for drugs that are delivered orally (for example, via a pill) because the active ingredient [needs to be absorbed into the lipids](https://emerypharma.com/blog/drug-lipophilicity-and-absorption-a-continuous-challenge-toward-the-goal-of-drug-discovery/) of biological membranes.

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



Which represents an 80:10:10 train:validate:test split


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
model.fit(train_dataset, nb_epoch=10)
```

    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_14:0", shape=(443,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_13:0", shape=(443, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Cast_4:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_17:0", shape=(2774,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_16:0", shape=(2774, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Cast_5:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_20:0", shape=(2385,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_19:0", shape=(2385, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Cast_6:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_23:0", shape=(188,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_22:0", shape=(188, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Cast_7:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_11:0", shape=(443,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_10:0", shape=(443, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_13:0", shape=(2774,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_12:0", shape=(2774, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_1:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_15:0", shape=(2385,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_14:0", shape=(2385, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_2:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_17:0", shape=(188,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_16:0", shape=(188, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_3:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_19:0", shape=(0,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_18:0", shape=(0, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_4:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_21:0", shape=(0,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_20:0", shape=(0, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_5:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_23:0", shape=(0,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_22:0", shape=(0, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_6:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_25:0", shape=(0,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_24:0", shape=(0, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_7:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_27:0", shape=(0,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_26:0", shape=(0, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_8:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_29:0", shape=(0,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_28:0", shape=(0, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_9:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_14:0", shape=(443,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_13:0", shape=(443, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Cast_4:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_17:0", shape=(2774,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_16:0", shape=(2774, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Cast_5:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_20:0", shape=(2385,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_19:0", shape=(2385, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Cast_6:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_23:0", shape=(188,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_22:0", shape=(188, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Cast_7:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_14:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_13:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Cast_4:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_17:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_16:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Cast_5:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_20:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_19:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Cast_6:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_23:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Reshape_22:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool_1/Cast_7:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_11:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_10:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_13:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_12:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_1:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_15:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_14:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_2:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_17:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Reshape_16:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_conv_1/Cast_3:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_14:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_13:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Cast_4:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_17:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_16:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Cast_5:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_20:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_19:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Cast_6:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(
    /usr/local/lib/python3.8/dist-packages/tensorflow/python/framework/indexed_slices.py:444: UserWarning: Converting sparse IndexedSlices(IndexedSlices(indices=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_23:0", shape=(None,), dtype=int32), values=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Reshape_22:0", shape=(None, 64), dtype=float32), dense_shape=Tensor("gradient_tape/private__graph_conv_keras_model/graph_pool/Cast_7:0", shape=(2,), dtype=int32))) to a dense Tensor of unknown shape. This may consume a large amount of memory.
      warnings.warn(





    0.6382771015167237



To check how well the model fits the train and test data, we examine the [Pearson correlation coefficient](https://www.scribbr.com/statistics/pearson-correlation-coefficient/) score.


```python
metric = dc.metrics.Metric(dc.metrics.pearson_r2_score)
print("Training set score:", model.evaluate(train_dataset, [metric], transformers))
print("Test set score:", model.evaluate(test_dataset, [metric], transformers))
```

    Training set score: {'pearson_r2_score': 0.35733596676636736}
    Test set score: {'pearson_r2_score': 0.15106334812602742}


As is typical, the score is better on the training set because the model fit itself to that data. By contrast, the test data is new to the model, so it hasn't had a chance to learn about that test data yet.

To get more information about how well the model works, let's fit it to the test data.


```python
lipos = model.predict_on_batch(test_dataset.X)
lipo_list = []
test_lipo_list = []
for molecule, lipo, test_lipo in zip(test_dataset.ids, lipos, test_dataset.y):
    lipo_list += [lipo[0]]
    test_lipo_list += [test_lipo[0]]
print(f"{lipo_list=}")
print(f"{test_lipo_list=}")
```

    lipo_list=[0.3845331, 0.3720408, -0.023484215, 0.531541, 0.77685773, 0.6073245, 0.2681268, 0.3807276, 0.26924998, 1.3168421, -1.1141771, 0.74688554, 1.2404268, 0.117568895, -0.2170724, 1.6985656, 1.1216733, 0.097339496, 0.8866433, 1.0619586, 0.023231909, 0.56239676, -0.6390465, 0.7374431, 0.621522, 0.50421333, 0.64951164, 0.58112335, 0.1416981, 0.94627386, 0.0364978, 0.8179041, -0.36597848, -0.3082065, 0.37811536, 0.9935026, 0.13907911, 0.06857388, -0.14017858, 0.23282218, 0.4666354, 0.91291934, 1.6561587, 0.8330316, -0.85378474, -0.0061901137, 0.9200858, 0.5301647, 0.14645682, 1.5272315, 0.6357204, 0.39473367, 1.2356297, 0.3559063, 0.5527899, 0.17531888, 0.928601, 1.1793733, 0.76913285, 0.397184, 1.0463703, 0.084844604, 0.8360811, 0.38981336, 1.1470077, 0.86523134, -0.06761201, 1.108468, 0.30433506, -0.7181876, 0.2707923, -0.1086918, 0.788031, 0.5218033, 1.0623503, 1.1002233, 0.72287184, 1.1593384, 0.89388514, 0.039959148, 1.0131426, 0.50972897, 0.8086295, 1.5717392, -0.64395314, 0.6700183, -0.1667041, -0.29344893, -0.014516905, 1.1411647, 0.9707738, 0.8250957, 0.7833464, 1.0754966, 0.91077024, -0.08869447, 1.6431936, 0.51833785, -0.20456116, 0.60434866, 0.41629457, 0.49047327, 0.294905, 1.1205404, 0.21038361, 1.0173624, 1.1023887, 1.0529292, 1.2233436, 0.8559557, 0.24101852, 0.41189367, 1.9521598, 0.1686906, -0.0029848963, 1.1588761, -0.093100294, -0.26053977, -0.8124754, 0.082770035, -0.08289644, 0.47479367, 1.2148663, 0.24129184, 1.2118862, 0.59788, 0.7221149, -0.85668874, -0.063501924, 0.40031344, 1.0817057, 0.8054461, 0.52570844, -0.39093864, 0.3379752, -0.036795944, 0.26691252, -0.47499305, 1.1067734, 0.48921388, 0.82112336, 0.7051856, -0.47432375, 0.62448794, 0.77601445, 0.790367, 0.5147836, 0.17237134, 0.05488892, 0.60998774, 0.19483237, 1.6400048, 0.54017365, 1.2681832, 0.53907734, 0.44269794, 0.40998566, 0.7260927, -0.92489886, 0.77017665, 0.90192366, -0.017482981, 0.5538161, 0.6853645, 0.95737594, -0.2036082, 0.5700011, 0.6709682, 1.0676379, 0.24085085, 1.4416652, 0.36321127, 0.97430956, 0.1699662, 0.9720545, -0.6201452, 0.12553944, -0.2703575, 0.8576147, 1.2413399, 0.78098714, 0.3072911, 1.4025961, 0.49760216, 0.8406737, 0.5681022, 0.87741125, 0.8626354, 0.48997545, 0.43210894, 0.40949142, 0.73812896, 0.9148929, -0.14690965, 0.26857364, 0.21870442, 0.3250382, 1.549915, -0.04997824, 1.0273519, 1.7752699, 0.5209258, 0.8910035, 0.6719488, 0.8051346, 0.05884485, 0.67517734, 0.7356062, 0.3399747, 0.41685158, 1.1449597, 1.2129986, 0.53483456, 1.0064216, 1.294542, 0.03149034, 0.031779364, 0.79296595, 1.1634603, -0.44887978, 0.4917606, -0.20404261, 1.5450358, 0.10812532, 1.32985, 0.8425587, 1.2763286, 0.13573556, 0.4076863, 0.2109316, 1.043422, 0.111568764, 0.48615992, 1.26736, 0.17908795, 0.5547333, 0.6477734, 0.71086156, 0.18833037, -0.05404, 0.1721112, 0.905608, 0.65089184, -0.3454305, 0.3881594, 0.18143158, 0.0060474128, 1.3791031, 1.0880893, -0.25525838, -0.38558608, 0.47776413, 0.5428214, -0.21291296, 0.5056131, 0.8428481, 0.24910347, 0.1737104, 0.51376575, 0.9410614, -0.27332973, 0.248318, 0.37836802, -0.0061704963, 0.80030054, 0.89736795, 1.3481767, 0.15114737, 0.65857816, 0.92081404, 0.557515, 0.34618044, 1.6231081, -0.51962894, -0.332995, 1.1102779, -0.29812247, 0.6740931, 1.2516999, -0.08390684, 0.90307033, -0.15123118, 1.1055475, 0.8578489, 1.1074977, 1.3402737, -0.10647899, 0.26910877, 0.67084694, 0.3697725, 0.118494496, 0.63540995, 0.576815, 0.7555401, 0.8233178, 0.4707631, -0.49241376, 0.16354932, 0.21740632, 0.8149489, 0.5347059, 0.9246714, 0.89023757, 0.2500642, 0.24600746, 1.0127708, -0.3351727, -0.4174571, 0.08440326, 0.3487339, 0.20730823, -0.08346535, 0.69527024, 0.7249452, 1.1218032, -0.47907853, -0.019051775, 0.56808025, 1.0376842, 0.7040532, 0.96236265, 0.82134223, 1.1633357, 0.81316423, 0.10947798, 0.19478966, 0.810393, -0.079371825, 0.47209954, 0.3798802, 0.89422375, 0.70015585, 1.0029, 0.2709477, -0.68354595, 1.2293369, 0.61893785, 1.0653851, 1.0299349, 1.5382998, 0.5256334, 1.092157, 0.4729516, 0.6358042, 0.40501332, 0.8597997, -0.49218243, 1.0113133, 0.7837634, 0.9936106, -1.0607456, 1.0815518, -0.078771606, 0.10964738, 0.39656943, -0.19871743, 0.90138286, 0.32232928, 0.4753825, 0.46027696, 0.22782572, 0.725273, 1.0067225, 1.0853257, -0.83870786, 0.7142333, -0.62045735, 0.027868733, 1.7499173, 0.8338685, 1.0742753, 0.21641599, 0.22085898, 0.8581599, 1.3379571, 0.4324636, -0.5490885, 1.435365, 0.91113853, 0.4992556, 0.7774311, 0.36286825, 0.7028681, -0.15022714, 0.61970806, 0.09741242, 0.48232627, -0.38753057, -0.18108341, 0.18954195, 1.501076, 0.32776356, 0.43698627, 1.1009026, 0.7166859, -0.27594167, 1.3813462, -0.5297604, 0.5254271, 0.92199063, 1.4334882, -0.31532514, 0.85784143, 0.277955, 0.062622264, -0.07856025, -0.06691642, 0.6150302, -0.09996809, 0.7752055, 0.8467773, -0.29572105, 0.38563007, 0.91909605, 0.6488579, 0.79765564, -0.22066335, 0.6076051, 0.1484554, 0.6295523]
    test_lipo_list=[-1.8108321895679695, 0.31965114479595336, -0.1923254704387878, 0.9389776954831404, 0.8564008220581824, 0.58389713975582, -0.6382405869335624, -0.6299828995910667, 0.01411671312360804, 0.6912470752082657, -0.3574792172887043, 0.29487808276846605, 1.5995926828828064, 0.6417009511532907, -0.4317984033711668, 1.3518620626079316, 0.6995047625507613, -0.6960443983310332, 1.3518620626079316, 1.054585318278082, -0.803394333783479, -0.2170985324662753, 1.8473233031576817, 0.195785834658516, -1.2575671376207493, -0.26664465652125013, 0.5426087030433406, -0.46482915274114994, -0.753848209728504, -0.10149090967133363, 0.44351645493339104, 1.43443893603289, -1.042867266715858, -1.951212874390399, -0.373994591973696, -0.2501292818362585, -0.5886944628785875, -2.85955848206494, 0.526093328358349, 0.03888977515109536, 1.4179235613478984, 0.8481431347156864, 0.8564008220581824, 1.3518620626079316, -0.6299828995910667, -0.13452165904131697, -0.48960221476863747, 0.5178356410158533, 0.44351645493339104, 0.8564008220581824, -0.2583869691787543, 0.13798202326104536, 1.104131442333057, -0.3822522793161918, -0.5143752767961249, -0.3822522793161918, 0.8151123853457031, -0.6382405869335624, -0.04368709827386298, 1.054585318278082, -0.018914036246375294, -0.3822522793161918, 0.575639452413324, -0.6299828995910667, 1.3518620626079316, 0.11320896123355768, 0.5426087030433406, 1.3188313132379483, -1.2080210135657745, 1.054585318278082, -0.3657369046312001, 1.43443893603289, -1.1171864527983204, 0.3526818941659367, 1.5748196208553193, 0.58389713975582, 1.1867083157580154, 0.8564008220581824, 0.06366283717858305, -1.1502172021683037, 0.36919726885092835, 0.30313577011096166, 1.43443893603289, 1.194966003100511, -0.8694558325234455, 0.34442420682344105, -1.5383285072656074, 1.244512127155486, -0.2996754058912335, -0.7125597730160248, -0.2170985324662753, 0.6086702017833073, 0.4270010802483994, -1.0758980160858413, 1.1619352537305276, -0.7951366464409831, 0.8564008220581824, -0.018914036246375294, 1.2692851891829737, 0.3113934574534577, 0.526093328358349, 1.0215545689080987, -1.232794075593262, 0.6086702017833073, -0.08497553498634197, 0.6417009511532907, -0.5308906514811166, 0.5178356410158533, -1.41446319712817, 0.4765472043033744, 0.44351645493339104, -0.05194478561635863, 1.1867083157580154, 1.0215545689080987, 1.2858005638679655, 0.6995047625507613, -2.2980357427752236, -0.0023986615613836274, 0.27836270808347435, 1.5500465588278316, -0.9768057679758914, 0.08017821186357435, 0.7242778245782491, 0.0719205245210787, 0.48480489164587004, 0.04714746249359138, -0.7951366464409831, -0.6960443983310332, 0.9389776954831404, -0.46482915274114994, -0.34096384260371265, 0.526093328358349, 0.022374400466103693, -1.4887823832106326, 1.4757273727453695, -0.2170985324662753, 0.6912470752082657, -0.803394333783479, 0.9967815068806114, 1.3518620626079316, -0.7951366464409831, 0.27836270808347435, -0.7868789590984873, 1.005039194223107, 0.21230120934350769, 0.7985970106607113, -0.0023986615613836274, -0.39876765400118347, -1.1171864527983204, -0.04368709827386298, 0.04714746249359138, 1.054585318278082, -0.7951366464409831, 0.31965114479595336, -0.13452165904131697, 0.58389713975582, 0.4022280182209117, 0.44351645493339104, -1.5465861946081032, 0.8151123853457031, 0.36919726885092835, 0.8233700726881991, -0.9107442692359247, -0.6299828995910667, 0.36919726885092835, 0.36919726885092835, 0.6004125144408117, 0.2618473333984827, -1.3484016983882035, 0.11320896123355768, 1.5170158094578483, 0.9224623207981487, 1.104131442333057, 1.3518620626079316, 1.2692851891829737, -1.951212874390399, -2.1163666212403154, 0.6086702017833073, 0.765566261290728, 0.6912470752082657, 0.4930625789883657, -1.224536388250766, 0.7077624498932573, 0.34442420682344105, 0.3774549561934244, -1.232794075593262, 1.1619352537305276, -0.852940457838454, 1.384892811977915, 0.3609395815084327, -0.3822522793161918, -0.14277934638381282, -0.5556637135086041, -0.10974859701382965, 0.6417009511532907, 1.4509543107178817, -0.4317984033711668, 1.43443893603289, 0.195785834658516, 0.9389776954831404, 1.682169556307765, -2.4466741149401487, -1.373174760415691, 1.104131442333057, 0.674731700523274, 1.054585318278082, 0.4022280182209117, -0.44005609071366264, 1.8142925537876984, 0.48480489164587004, 0.44351645493339104, 1.6904272436502608, -0.5061175894536292, 1.384892811977915, 0.9389776954831404, -0.1923254704387878, -0.7786212717559915, 0.06366283717858305, 1.5995926828828064, -0.2996754058912335, -2.0172743731303653, 1.0215545689080987, 0.6086702017833073, 0.6995047625507613, 1.9299001765826398, -0.2170985324662753, 1.682169556307765, -0.5556637135086041, -0.1180062843563253, -0.7951366464409831, 0.8233700726881991, -1.1006710781133286, 0.6995047625507613, 0.41874339290590334, -0.13452165904131697, 0.773823948633224, 0.6086702017833073, 1.8968694272126563, -0.10974859701382965, -1.0263518920308663, -0.423540716028671, 0.46003182961838274, -0.05194478561635863, -1.5300708199231117, 0.20404352200101203, -1.1254441401408162, -0.2170985324662753, 1.682169556307765, -0.02717172358887131, -1.3979478224431785, 0.8151123853457031, 0.01411671312360804, 1.0876160676480653, 0.7077624498932573, 1.3518620626079316, 1.294058251210461, -0.22535621980877116, -0.7951366464409831, -1.1171864527983204, 0.3609395815084327, 1.43443893603289, 0.773823948633224, -0.20058315778128363, -1.183247951538287, -0.31619078057622513, 0.44351645493339104, 0.8564008220581824, 1.0958737549905613, -0.04368709827386298, 0.6499586384957867, 1.0215545689080987, 0.05540514983608703, 0.1627550852885327, 0.5508663903858367, -1.951212874390399, 1.7069426183352527, -0.6464982742760583, 1.2858005638679655, 0.526093328358349, 0.12972433591854934, 1.104131442333057, -0.8281673958109664, 0.44351645493339104, -1.0924133907708329, -0.3574792172887043, 0.8233700726881991, -0.993321142660883, -1.034609579373362, 0.5426087030433406, 0.6912470752082657, 1.5170158094578483, -0.7290751477010164, -1.323628636360716, 1.0711006929630738, 0.6664740131807784, 0.955493070168132, 0.31965114479595336, 0.7325355119207447, 0.9389776954831404, 1.2610275018404777, 0.9389776954831404, -0.0023986615613836274, 0.27836270808347435, 0.195785834658516, -0.5804367755360916, 0.4765472043033744, -1.1419595148258077, -2.132881995925307, 1.145419879045536, 0.2205588966860037, -0.9850634553183871, -1.1667325768532952, 0.8564008220581824, -0.3822522793161918, -0.184067783096292, 0.7077624498932573, -2.5292509883651064, -1.2080210135657745, 0.7985970106607113, 1.1619352537305276, 0.1214666485760537, 0.9389776954831404, 0.15449739794603704, 0.9389776954831404, -1.1089287654558246, 0.8564008220581824, -0.2996754058912335, -1.4970400705531284, -0.852940457838454, 0.30313577011096166, -0.15103703372630864, -0.8777135198659413, 0.1627550852885327, 1.2032236904430071, -0.9520327059484038, 1.013296881565603, 0.09669358654856601, 1.4261812486903944, -0.423540716028671, 0.8481431347156864, -0.5556637135086041, -0.5474060261661083, 0.03063208780859971, 1.624365744910294, 0.1214666485760537, 0.28662039542597, -1.3318863237032117, 0.8564008220581824, 1.1206468170180488, 0.5426087030433406, -0.5061175894536292, 0.6912470752082657, -0.7125597730160248, -0.7455905223860081, -0.4565714653986541, -1.1749902641957912, 0.7820816359757197, -0.3822522793161918, 0.4270010802483994, 0.3774549561934244, -0.48134452742614164, 0.7490508866057364, 0.7490508866057364, -0.9602903932908996, -1.2080210135657745, -0.13452165904131697, 0.8564008220581824, -1.133701827483312, 0.963750757510628, 0.195785834658516, 0.7573085739482324, 1.1289045043605443, -2.2484896187202486, 0.8564008220581824, 1.9133848018976485, 0.27836270808347435, -0.8694558325234455, 0.7903393233182157, 1.8473233031576817, 0.765566261290728, 1.7482310550477318, -1.6209053806905658, -0.10149090967133363, 0.03063208780859971, 1.2692851891829737, -0.423540716028671, 1.2858005638679655, -1.0511249540583536, 0.3609395815084327, -0.13452165904131697, 0.25358964605598705, -0.6052098375635792, -0.5061175894536292, 0.6912470752082657, 0.7985970106607113, 1.0215545689080987, 0.58389713975582, 0.44351645493339104, -0.803394333783479, -0.373994591973696, 1.3518620626079316, -0.13452165904131697, 1.0793583803055697, 0.09669358654856601, -0.8777135198659413, 0.03063208780859971, 0.6169278891258033, -0.6299828995910667, 0.31965114479595336, 1.7895194917602106, -0.08497553498634197, 0.3609395815084327, -0.8446827704959582, 0.08843589920607037, 1.4592119980603777, 0.8151123853457031, 1.3270890005804443, -0.17581009575379614, 0.0719205245210787, 0.8068546980032074]


Then we make put the measured and model-predicted results into a [pandas dataframe](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) for easy processing.


```python
import pandas as pd
df = pd.DataFrame(list(zip(test_lipo_list, lipo_list)), columns = ["measured", "predicted"])
df
```





  <div id="df-0e7841ee-894c-419c-9c0d-71890a18a6a1">
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
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>-1.810832</td>
      <td>0.384533</td>
    </tr>
    <tr>
      <th>1</th>
      <td>0.319651</td>
      <td>0.372041</td>
    </tr>
    <tr>
      <th>2</th>
      <td>-0.192325</td>
      <td>-0.023484</td>
    </tr>
    <tr>
      <th>3</th>
      <td>0.938978</td>
      <td>0.531541</td>
    </tr>
    <tr>
      <th>4</th>
      <td>0.856401</td>
      <td>0.776858</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>415</th>
      <td>0.815112</td>
      <td>0.797656</td>
    </tr>
    <tr>
      <th>416</th>
      <td>1.327089</td>
      <td>-0.220663</td>
    </tr>
    <tr>
      <th>417</th>
      <td>-0.175810</td>
      <td>0.607605</td>
    </tr>
    <tr>
      <th>418</th>
      <td>0.071921</td>
      <td>0.148455</td>
    </tr>
    <tr>
      <th>419</th>
      <td>0.806855</td>
      <td>0.629552</td>
    </tr>
  </tbody>
</table>
<p>420 rows × 2 columns</p>
</div>
      <button class="colab-df-convert" onclick="convertToInteractive('df-0e7841ee-894c-419c-9c0d-71890a18a6a1')"
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
          document.querySelector('#df-0e7841ee-894c-419c-9c0d-71890a18a6a1 button.colab-df-convert');
        buttonEl.style.display =
          google.colab.kernel.accessAllowed ? 'block' : 'none';

        async function convertToInteractive(key) {
          const element = document.querySelector('#df-0e7841ee-894c-419c-9c0d-71890a18a6a1');
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




Now we can use a scatter plot to compare the predicted against measured values. We use the [seaborn statistical data visualization package](https://seaborn.pydata.org/) to easily plot the data.


```python
import numpy as np
[np.min(lipo_list), np.max(lipo_list)]
```




    [-1.1141771, 1.9521598]




```python
_[1] - _[0]
```




    3.0663369




```python
[np.min(test_lipo_list), np.max(test_lipo_list)]
```




    [-2.85955848206494, 1.9299001765826398]




```python
_[1] - _[0]
```




    4.78945865864758




```python
import seaborn
```


```python
seaborn.scatterplot(data=df, x = "measured", y = "predicted").set(title='Lipophilicty for test data:\noctanol/water distribution coefficient\n(logD at pH 7.4)\n')
```




    [Text(0.5, 1.0, 'Lipophilicty for test data:\noctanol/water distribution coefficient\n(logD at pH 7.4)\n')]




    
![png](2022-12-12-DeepChem_lipophilicty_machine_learning_files/2022-12-12-DeepChem_lipophilicty_machine_learning_25_1.png)
    


One observation is that the the measured lipophilicty data have a greater range (about -2.9 to 1.9, comprising 4.8 units) than the predicted results (about -1.1 to 2.0, comprising 3.1 units) by about 1.7 units. This suggests that outliers in the lipophilicty data are not well modeled.

Let's compare to the train data set.


```python
lipos_train = model.predict_on_batch(train_dataset.X)
lipo_list_train = []
test_lipo_list_train = []
for molecule, lipo, test_lipo in zip(train_dataset.ids, lipos_train, train_dataset.y):
    lipo_list_train += [lipo[0]]
    test_lipo_list_train += [test_lipo[0]]
df_train = pd.DataFrame(list(zip(test_lipo_list_train, lipo_list_train)), columns = ["measured", "predicted"])
seaborn.scatterplot(data=df_train, x = "measured", y = "predicted").set(title='Lipophilicty for train data:\noctanol/water distribution coefficient\n(logD at pH 7.4)\n')
```




    [Text(0.5, 1.0, 'Lipophilicty for train data:\noctanol/water distribution coefficient\n(logD at pH 7.4)\n')]




    
![png](2022-12-12-DeepChem_lipophilicty_machine_learning_files/2022-12-12-DeepChem_lipophilicty_machine_learning_28_1.png)
    



```python
import numpy as np
[np.min(lipo_list_train), np.max(lipo_list_train)]
```




    [-2.2968836, 2.1289487]




```python
_[1] - _[0]
```




    4.78945865864758




```python
[np.min(test_lipo_list_train), np.max(test_lipo_list_train)]
```




    [-3.0247122289148565, 1.9299001765826398]




```python
_[1] - _[0]
```




    4.954612405497496



The train data also has a wider measured range (about -3.0 to 1.9, comprising 5.0 units) than predicted (about -2.3 to 2.1, comprising 4.8 units) by about 0.2 units, though the difference in ranges is quite a bit less than the 1.7 units differential for the test data.

While the range differential seems like a minor factor in the train set, we might want to change how we featurize the data to capture a wider range of molecular characteristics. We used the [`GraphConv`](https://deepchem.readthedocs.io/en/latest/api_reference/featurizers.html#graph-convolution-featurizers) featurizer, which represents only the *atoms* in a molecule. We might consider the [`WeaveFeaturizer`](https://deepchem.readthedocs.io/en/latest/api_reference/featurizers.html#weavefeaturizer) which also represents the bonds, though it requires more resources.

Perhaps a bigger factor is as follows. We note that both the [Pearson correlation coefficient](https://www.scribbr.com/statistics/pearson-correlation-coefficient/) score is worse, and the range differential is so much greater, for the test data than for the train data. One reason might be that our train set might not have molecules similar to those in the test set. Recall that we split by scaffold, so it's possible that splitting led to compounds in the test set that have scaffolds significantly different from those in the train set. Adding compounds to the dataset so that there is less "scaffold distance" (difference in scaffold structure) between groups of compounds might help.

To take a simple hypothetical example, if we had a dataset with compounds containing either one or three rings, the scaffold splitting might put all the one-ring compounds in the train set and all the three-ring compounds in the test set. We expect that would lead to poor predictions on the test set. We would want to augment the dataset by adding compounds containing two rings.

This blog post was based on the DeepChem tutorial [The Basic Tools of the Deep Life Sciences](https://github.com/deepchem/deepchem/blob/master/examples/tutorials/The_Basic_Tools_of_the_Deep_Life_Sciences.ipynb).

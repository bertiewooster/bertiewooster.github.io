#!/bin/bash

# Note: Before running this script,
#  1) If you created the notebook locally
#     a) format with Ruff by ensuring no cell is selected, then pressing Option+Shift+F on Mac (Alt+Shift+F on Windows)
#     b) You can make it into a Google Colab by
#        i) creating a new Google Colab on bertiewoostersharing@gmail.com account
#        ii) choosing File > Upload notebook
#  2) If you downloaded the notebook from Google Colab, replace file name underscores with hyphens
#  3) If you want to embed an image in the notebook rather than importing it from an external file,
#     a) Set up Markdown cell with <img src="" />
#     b) Convert image to Base64 at https://www.base64-image.de/ (can be JPG, PNG, GIF, WebP, SVG or BMP, up to 1MB)
#     c) Paste the <img> output (starting with data:image/) into the quotes in the img tag
#  4) Change the name of the notebook (input) and markdown (output) files in the non-commented code below
#  5) Save this script
#  6) Use venv in this folder, not conda, in terminal using 
#        . venv/bin/activate

# Run this script in the terminal using:
# ./nb_md.sh

jupyter nbconvert _notebooks/Cheminformatics-Blogging-How-To.ipynb --to markdown --output Cheminformatics-Blogging-How-To.md

# Note: After running this script, to post a blog post:
#  1) Move the markdown file into the _posts folder
#  2) Move the _files folder (containing the images) into the images folder
#  ----- The rest of these instructions are for the markdown file -----
#  3) Fix the image links in the markdown file by converting e.g.
#          ![png](2022-11-11-RDKit-Recap-decomposition-tree_files/2022-11-11-RDKit-Recap-decomposition-tree_10_0.png)
#      to
#         ![png](/images/2022-11-11-RDKit-Recap-decomposition-tree_files/2022-10-09-RDKit-find-and-highlight-the-maximum-common-substructure-between-molecules_24_0.png)
#      by replacing
#         ![png](
#      with
#         ![png](/images/
#      a) If there are any svg images, do likewise-- replace
#             ![svg](
#         with 
#            ![svg](/images/
#  4) If you have a link to an image in the top-level /images folder, in the markdown file, replace
#            ../images
#         with 
#            /images
#  5) Add alt tags to Jupyter images by replacing ![png] or ![svg] with a description, e.g. ![Annotated tree]
#  6) If notebook had molframe graphs, delete from
#            <style>*[data-root-id]
#         to the </div> after
#            float: left'></img>
#
#         and from
#            <div id='p
#         to
#            })(window);</script>
#         (You can triple-click on the enormous block of text to select it)
#  ---If you downloaded the notebook from Google Colab; if not, create the Google Colab version and do the opposite to it---
#  7) Delete the first paragraph:
#            *Blog post by [Jeremy Monat](https://bertiewooster.github.io/)*
#  8) Add, to the end of the first section:
#            *[Open this notebook in Google Colab](link) so you can run it without installing anything on your computer*


# Rating Survey Beginning to End, Part 3: Publish the results online

[Part 1: Create a survey in Google Forms, collect responses, then export the responses to a .csv file](/2022/03/20/survey-pipeline-1.html)

Part 2: Import the survey response data into a GitHub repository, then analyze the data

[Part 3: Publish the results online using GitHub Pages](/2022/03/26/survey-pipeline-3.html)

## Part 2: Analyzing data

## Forking and cloning the GitHub repository
You can get all the code needed for the analysis:
1. [Fork the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) at https://github.com/bertiewooster/survey-analysis.
2. The forked repository will be *username*/survey-analysis where *username* is your GitHub username.
3. [Clone the forked repository](https://docs.github.com/en/enterprise-cloud@latest/repositories/creating-and-managing-repositories/cloning-a-repository) onto your computer. By default, it will be cloned onto your computer into a folder called survey-analysis.

## Copy in the survey data
Copy the rating-options-survey.csv file into the survey-analysis folder, replacing the existing rating-options-survey.csv file.

## Analyze the data
TODO How to customize Jupyter Notebook for their data

## Export the notebook to HTML
1. In your terminal, go to the survey-analysis folder. If you're using [Visual Studio Code, pull up the terminal](https://code.visualstudio.com/docs/editor/integrated-terminal) and it will be in that folder.
2. Run `./nb_html`
3. Two HTML versions of the notebook will be created in the doc folder:
    1. Analysis.html has only the results from the notebook. This is the cleaner version that I recommend you publish.
    2. Analysis-with-code.html: If you want to show your community the code used to analyze the results, you can publish this longer version.

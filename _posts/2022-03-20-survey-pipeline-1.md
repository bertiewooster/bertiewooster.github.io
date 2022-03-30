# Rating Survey Beginning to End, Part 1: Create a survey, collect and export data

## Free tools for survey creation, analysis, and publishing of results

This tutorial is at the intermediate level: each step that's specific to this workflow is explained, and steps that are commonly used elsewhere will have instructions linked to on other sites.

Part 1: Create a survey in Google Forms, collect responses, then export the responses to a .csv file

[Part 2: Import the survey response data into a GitHub repository, then analyze the data](/2022/03/26/survey-pipeline-2.html).

[Part 3: Publish the results online using GitHub Pages](/2022/03/26/survey-pipeline-3.html)

## Part 1: Creating the survey, collecting and exporting data

## Scenario

To make a choice based on community input, you may need to ask a group of people to rate a set of options.

## Why use Google Forms, Jupyter Notebook, and GitHub Pages
- All of them are free
- The output from one tool can be easily used in the next tool in the chain
- Google Forms is an easy, free way to create a survey
- Google Forms can easily export survey data
- Jupyter notebooks are easy to use and inspect; host narratives, code, and outputs; and can produce high-quality graphs
- Jupyter notebooks can be easily converted to HTML
- GitHub Pages is an easy way to post HTML files

## What the final result will include
After you follow these instructions, you'll have a web page (and Jupyter Notebook) that includes:
- A table of the mean and standard deviation of the rating for each option
- A bar graph of the mean and standard deviation of the rating for each option
- A histogram for each option with the count of responses for each rating level, including a dashed vertical line at the mean rating

TODO: Include those results from survey-analysis repo when done

## What you'll need
- A Google account
    - You can [create a Google account for free](https://accounts.google.com/signup)
- The ability to edit Jupyter Notebooks
    - You can [follow this tutorial to set up a Jupyter Notebook](https://www.dataquest.io/blog/jupyter-notebook-tutorial/)
    - Or you can edit [Jupyter Notebooks in Visual Studio Code](https://code.visualstudio.com/docs/datascience/jupyter-notebooks).
- A GitHub account
    - You can [create a GitHub account for free](https://github.com/join)

## Setting up a survey in Google Forms
When asking for feedback on a set of options, two options are rating and ranking. Rating means giving each option a score on a scale, for example "How helpful is this option" on a 1-4 scale from "Not very helpful" to "Very helpful". Ranking means putting the options in order of preference.

Asking respondents to rate a set of options is a good choice:
- Has a relatively low cognitive load because the respondent need keep only one option in mind at a time.
- To rank the options, the respondent would have to have all the options in mind at once, then decide which option is best, second best, etc.
- Forcing respondents to order the options may also distort their feedback on how good the options are. For example, a respondent might think options A and B are nearly as good, yet be forced to rate one as best and the other as second best.

The structure of the rating survey is simple:
- Introductory page: Why this feedback is requested
- One page for each option
    - Rate the option on the scale
    - What do you LIKE about the option?
    - What do you DISLIKE about the option?
- Final page: Do you have any other comments?

Clone, edit, and share the survey by:
1. [Clicking this link to open the provided Rating options survey folder](https://docs.google.com/forms/d/1Z7atj9bw9FsVih4bs_XgUmlHpsbYBO3UAOoUuKzoUIo/edit?usp=sharing) in your Google Drive
    - If you haven't already logged in to your Google account, you will be prompted to log in
2. Right-clicking on the Rating options survey form and selecting Make a copy
<img src="/images/copy-google-form.png" style="border: 1px solid gray; width: 50%; display:block;">
3. Clicking Show file location
<img src="/images/google-drive-show-file-location.png" style="border: 1px solid gray; width: 50%; display:block;">
4. Renaming from "Copy of Rating options survey" to the name you'd like for your survey, for example "Rating options survey"
5. Double-clicking the form you just renamed
6. [Editing the Google form](https://support.google.com/docs/answer/2839737?hl=en) to ask the questions you need answered
7. Sharing your form using the Send button at the top right of your form
    - You can share by email, getting and posting a link, embedding HTML, or via Facebook or Twitter

A week or two is typically a good length of time to gather responses.

## Exporting survey data from Google Forms
1. Click Responses at the top of the form
<img src="/images/google-form-responses-3dots.png" style="border: 1px solid gray; width: 50%; display:block;">

2. Click the three-dot icon

3. Choose Download responses (.csv)
<img src="/images/google-form-download-responses.png" style="border: 1px solid gray; width: 50%; display:block;">

4. Copy survey-responses.csv into the GitHub repo.

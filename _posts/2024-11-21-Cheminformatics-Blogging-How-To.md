# How to Write Cheminformatics Blog Posts

As the YouTubers would say, “A lot of you have been asking me about how to write cheminformatics blog posts.” Well, not a lot, but at least a couple! I realized that there’s a pattern to how I write cheminformatics blog posts (16 so far), so I’m sharing that here.

My blog posts are intended to be tutorials that explore a topic, usually with existing tools. I figure out how to accomplish a cheminformatics objective, then share that.

This post is about the creative process more than the technical details and is pitched toward people with a chemistry background.

Throughout this post, I’ll use as an example my last blog post, [Why some organic molecules have a color: Correlating optical absorption wavelength with conjugated bond chain length]]({% post_url 2024-10-15-Color-from-Conjugation %}).

## Find a jumping-off point

A jumping-off point can be any of these approaches you find interesting: 
- A tutorial. Examples: [Getting Started with the RDKit in Python](https://www.rdkit.org/docs/GettingStartedInPython.html), [RDKit Cookbook](https://rdkit.org/docs/Cookbook.html), or a cheminformatics blog post by, for example, [Takayuki Serizawa](https://iwatobipen.wordpress.com/).
- A talk or presentation. Examples: [RDKit User Group Meeting (UGM) recordings](https://www.youtube.com/playlist?list=PLugOo5eIVY3EHeBuSABISVok5-Q7kE0O1).
- A dataset. Examples: [Experimental database of optical properties of organic compounds](https://www.nature.com/articles/s41597-020-00634-8), datasets on [MoleculeNet](https://deepchem.readthedocs.io/en/latest/api_reference/moleculenet.html).

**Example:** My background is in optical spectroscopy, so I searched for an optical spectroscopy dataset, and found one with thousands of chromophores and their absorption and emission wavelength maxima. Thanks to relatively new publications such as [Nature *Scientific Data*](https://www.nature.com/sdata/), there are a good number of datasets available.

## Ask yourself "What would I as a chemist want to do with this?"

This is where your chemistry background comes into play. Given the jumping-off point, think about what you as a chemist would want to do with the function or topic.

Your interests will determine what you find interesting as a jumping-off point and what you want to do with it. It’s that combination that means that each cheminformatician will come up with different topics.

Find something compelling to you and follow your curiosity. When you find something you want to work on and can’t put down, you’ve found a good topic.

**Example:** When I considered the optical dataset, I remembered that there's a correlation between conjugated bond chain length and absorption wavelength. So I thought about verifying that using this dataset.

## Formulate an approach

### What chemical information do you want?

I think it’s best to start with the question of what chemical information you want, then figure out a way to calculate it. This order ensures that the topic is chemically relevant.

**Example:** I needed to determine the longest conjugated bond chain length in a molecule. (Absorption wavelength was already provided, though I wanted to convert it to energy, so that required a calculation.)

### What code will you need to write to extract the chemical information?

This is where the cheminformatics part comes in: Translating the chemical information question into a plan for how to get it using code. Sometimes the core function you need is already available, for example in a cheminformatics toolkit, so it's a matter of applying it to your situation. Other times, you need to write the function yourself, for example by considering a molecule as a graph where the nodes are atoms and the edges are bonds.

**Example:** The optical dataset provided SMILES, which the RDKit can turn into a molecular graph, so the task was to devise an algorithm that would traverse the conjugated bonds in the molecular graph and determine which were connected to a starting bond. Then find all such chains and keep the longest one in the molecule.

## Consult a journal article

It’s often helpful to consult a journal article to get molecular structures, experimental details, etc. I’ve often used my American Chemical Society membership benefit [ACS Member Universal Access (50 article accesses per term)](https://www.acs.org/membership/member-benefits.html) to access journal articles, though many recent publications may be available on [preprint archive servers](https://chemrxiv.org/engage/chemrxiv/public-dashboard) now.

**Example:** I found the optical properties dataset via the [journal article that presented it](https://www.nature.com/articles/s41597-020-00634-8). I also read the [journal article that studied the oligomers](https://pubs.acs.org/doi/full/10.1021/acs.joc.8b00311) in the homologous series. This helped me correct the structures in the database and understand the authors’ explanations of trends within that series.

## Code your way through the project

Here's where I start writing code. It's typically an iterative process where I solve one problem and get results, only to realize I need to solve an additional problem.

**Example 1:** In finding the longest conjugated chain, I started with one algorithm, but realized it didn’t produce the results that made chemical sense because it wouldn’t include the shared bond between the two rings in naphthalene. So I generalized my conception of a “chain” to include all conjugated bonds adjacent to an existing conjugated bond in the “chain,” which allowed for multiple paths.

**Example 2:** After that code was working, the results showed that there wasn't much of a correlation of absorption wavelength against longest conjugated bond chain length if I included thousands of chromophores. So I realized I needed to narrow down to a series of molecules with similar structures by finding chromophores in the dataset with moieties in common.

#### Make the critical code more efficient

Often, there's a single cheminformatics operation which takes a significant amount of time on a sizeable dataset. This slow step is worth speeding up, especially because I often run it dozens of times while drafting the blog post.

**Example:** In the optical properties post, the slow step was traversing each molecular graph to determine its longest conjugated bond chain. I sped this up by only traversing each bond once, rather than using each bond as a starting point and then traversing all connected bonds. I also sped up the operation at the dataframe level: Because it contained multiple rows for many chromophores, I cached the result for each chromophore so it wouldn't need to be recalculated.

Note that speeding up the code is not the only goal: Because the blog post is intended to be educational, there's a balance between code speed and clarity. For example, in many posts, I use Polars in eager mode rather than [lazy mode](https://docs.pola.rs/user-guide/lazy/) because I want to show the steps along the way, for example the dataframe after the data is initially loaded. If speed was the sole goal, I'd use lazy mode to evaluate only once, to get the final results.

#### Clean up the code

In getting the code working, I often go down paths that don't work out, and include debugging code. So I remove that code before posting. I also rename functions whose purpose has changed since I originally wrote them, and add docstrings.

**Example:** During development, I displayed the Polars dataframe to check the results of individual steps. Once those steps were working, I removed that debugging code.

## Write a clear narrative

It's easiest for the reader to understand the post if there's a clear narrative. So after the code is working, I start writing the narrative. I might start writing a narrative (Markdown text) block, then realize that the narrative would make more sense if I rearranged the code blocks, so I do that and then adjust the narrative.

**Example:** To find a homologous series, I used four substructure matches. During development, I ran all four at once. But to make the process clearer for the reader, I broke it out into three steps (with one, one, and two substructures, respectively).

## Review

Before publishing, clear out the existing results by restarting the Jupyter kernel or closing and reopening the notebook. Then re-run all the code in order to make sure there are no problems when it runs in top-to-bottom order.

Re-read the narrative along with the updated code outputs to make sure the descriptions are still correct.

## Publish!

There comes a time when the post is done (enough) and ready to publish. After publishing to a web page, I then review the published page and correct any broken images or links. Finally, I publicize the post on social media such as LinkedIn and Mastodon.

## Make any changes after publication

Readers may point out improvements you could make, or the data sources could be updated. Ideally, these can help you accomplish the same results with code that's shorter or more straightforward.

**Example:** After I mentioned the discrepancy in a few molecular structures to the maintainers of the optical properties dataset, they updated it. I pulled in the new dataset, which let me simplify the notebook by removing the manual correction of structures.

# Murcko Scaffolds Post on the RDKit Blog: An LLM Circle

## An LLM Chat Plants the Seed

To come up with small projects to do in preparation for job interviews (contact me if you need a cheminformatics scientific software developer!), I interactively asked ChatGPT for suggestions. It suggested

> **Cheminformatics-Flavored Algorithm Problem**
> A classic coding challenge but wrapped in a domain-relevant context.
>
> Example Prompt:
>
> "You are given a list of SMILES strings. Write a function that identifies and counts the number of unique scaffolds using the Bemis-Murcko scaffold definition. Use RDKit."

## The Murcko Scaffolds Blog Post

I wrote code to find the Murcko scaffold for a set of molecules, then determine the count of unique scaffolds. That answered the ChapGPT challenge.

But as usual, I wanted to visualize the results. So I grouped the molecules by scaffold, then used [`MolsMatrixToGridImage`](https://greglandrum.github.io/rdkit-blog/posts/2023-10-25-molsmatrixtogridimage.html) (which I contributed to the RDKit) to display each scaffold as a row with its molecules to the right.

I contributed the [Murcko Scaffolds Tutorial](https://greglandrum.github.io/rdkit-blog/posts/2025-06-06-murcko-scaffolds.html) post to the RDKit blog. Thanks to Greg Landrum for publishing the post on that blog!

![Molecular grid displaying each of three Murcko scaffolds as a row, with 1-3 molecules to the right of that scaffold](/images/murcko_scaffold_matrix.png)

## Another LLM Tool is Applied to Reproduce My Visualization

A couple weeks later, I saw a [LinkedIn post about an RDKit model context protocol (MCP) server](https://www.linkedin.com/feed/update/urn:li:activity:7343304678215970816/) whose molecular visualization looked awfully familiar: A row for each scaffold with its molecules to the right. The authors indeed used my blog post as a test case of their RDKit Copilot!

## Put Things Out Into the World, No Matter How Simple

I love that my simple code and visualization was an inspiration for others to test their new tool. It shows that you never know what effect it will have when you put something out into the world: Someone else may take it in a direction you never expected.

It also confirms what I've found in other ways: Simple ideas make good blog posts. It's easy to think that you have to do something very complex to make a contribution: People have asked me if they need to come up with very advanced ideas for blog posts. But as I've written in [How to Write Cheminformatics Blog Posts]({% post_url 2024-11-21-Cheminformatics-Blogging-How-To %}), coming up with ideas can be as simple as going through a tutorial and building on it. In fact, actually I've found via Google Search Console that my very first blog post, [Find and Highlight the Maximum Common Substructure Between a Set of Molecules Using RDKit] is the most visited from searches. Presumably this is because many people want to accomplish the simple task of finding a maximum common substructure.

So no matter how simple the idea, if you learn something by coding it, others probably will too, so go ahead and blog about it.

If your inspiration leads you to blog about something involved like [tautomers]({% post_url 2024-05-01-Tautomer-Sources-Comparison %}), that's fine too! Just don't think you should reject an idea because it's too simple.

## About

We provide examples of computing operations on genomic ranges using the
[tidy data philosophy](https://tidyr.tidyverse.org/articles/tidy-data.html).

You can find the preview of this example at 
<https://tidyomics.github.io/tidy-ranges-tutorial/>

---

## Contributing

We use quarto to publish the website. Computation is done locally through `quarto render` with quarto's `freeze: auto` option (see `_quarto.yml`) and the `_freeze/` directory is commited. The website is rendered through a github workflow (see `.github/workflows/publish.yml`)

If you want to contribute: 

1. fork this repository and create a feature branch. 
2. render the website using `quarto render` command to inspect the results locally (you need to serve `_book/`). 
3. Push feature branch changes to your github fork. The github workflow should work on your fork and feature branch, so one can inspect the rendered version hosted on your fork using github pages. If it does not work, make sure to go on your github forked repository and enable workflows.
4. Open a PR on main repository and share the url to the rendered website on your fork so we can inspect the rendered website.


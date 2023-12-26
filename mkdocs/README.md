
## Overview of Project Components

This file provides an overview of the different components in the project directory. It includes information about the following files and folders:

- `mkdocs.yml`: This file is the main configuration file for the MkDocs documentation project. It contains settings such as the site name, theme, and navigation structure.

- `package.json`: This file is used to save the bash commands for developing and building the project.

- `index.md`: This file serves as the main documentation page for the project. It provides an introduction and overview of the project's purpose and contents. It is analogus to the index.html file for websites.

- `img` folder: This folder contains images used in the documentation. It includes images that are referenced in the Markdown files to enhance the content. The `img` folder is also used for the site logo and favicon.


# How to develop
```bash
mkdocs serve
```

## Convert Jupyter to Markdown
You can use the nbconvert command-line tool to convert Jupyter notebooks to markdown. To do this, open a terminal window and navigate to the directory where the Jupyter notebook is located. Then, run the following command:
```bash
jupyter nbconvert --to markdown <jupyter_notebook_name>.ipynb
```
To convert all markdown files in `src` folder and put them in a output folder

```bash
find src -name "*.ipynb" -exec sh -c 'jupyter nbconvert --to markdown --output-dir=mkdocs/docs/2_notebooks/$(basename "$(dirname "$(dirname "{}")")")/$(basename "$(dirname "{}")") "{}"' \;
```
In this command, `$(basename "$(dirname "$(dirname "{}")")")` is used to extract the name of the second parent folder, and `$(basename "$(dirname "{}")")` is used to extract the name of the immediate parent folder. The converted Markdown files will be placed in the `mkdocs/docs/notebooks/` directory with both parent folder names included in the file name.

## Extract the key value pairs from the generated md files.

The files will be sorted first by the parent's parent directory and then by the parent directory within each group of the first sort. This nested sorting will reflect a two-level directory hierarchy in the files.txt.

```bash
find mkdocs/docs/2_notebooks -type f -name "*.md" -exec bash -c 'f="{}"; 
parent=$(dirname "$f"); 
grandparent=$(dirname "$parent"); 
basename=$(basename "$f"); 
echo "$(basename "$grandparent")/$(basename "$parent") - - $basename: ${f#mkdocs/docs/}"' \; | sort -t '/' -k1,1 -k2,2 | sed 's/^[^ ]*\/[^ ]* - - /- /' > mkdocs/files.txt
```




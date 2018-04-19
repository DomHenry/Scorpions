
# Notes on RMarkdown stuff for scorpion document --------------------------

# Notes
# 1. adding a figure caption in the chunk header moves figures to the top of the page
# (would need to follow the link below to fix)
# link:https://tex.stackexchange.com/questions/101725/latex-figures-appear-before-text-in-pandoc-markdown
# 
# 2. Options for ploting size
# https://sebastiansauer.github.io/figure_sizing_knitr/
#   https://www.zevross.com/blog/2017/06/19/tips-and-tricks-for-working-with-images-and-figures-in-r-markdown-documents/
#   http://blog.revolutionanalytics.com/2017/06/rmarkdown-tricks.html

# in kable() use argument 'escape = FALSE' to render maths equations in a table

# backslash followed by two whitespaces creates new line

# grid.arrange creates a blank page after plots. See error here:
# https://github.com/wilkelab/cowplot/issues/24

# One solution would be to use pdftk at the end to remove the page (but it will be 
# numbered!!! AAARRR)

#https://stackoverflow.com/questions/12481267/in-r-how-to-prevent-blank-page-in-pdf-when-using-gridbase-to-embed-subplot-insi


# The error to google here is" rmarkdown is not plotting within margins of pdf using loop"
# https://stackoverflow.com/questions/37030219/r-markdown-plots-within-a-loop-going-out-of-margin-when-typesetting-to-pdf

# See this page for more table options:
#   https://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf
# Comment out text shortcut is some in normal R

# https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html

# TEMP REMOVE TO IMPROVE RENDER TIME --------------------------------------


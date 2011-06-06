#!/bin/bash
# Export from Lyx to Latex
# Compile .tex to .pdf
if [ -f /Applications/LyX.app/Contents/MacOS/lyx ];
then
    lyx_executable="/Applications/LyX.app/Contents/MacOS/lyx"
else
    lyx_executable="lyx"
fi
echo $lyx_executable

$lyx_executable --export latex tutorial
pdflatex -shell-escape tutorial.tex
bibtex tutorial
pdflatex -shell-escape tutorial.tex

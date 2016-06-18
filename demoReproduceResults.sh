 #!/bin/sh

TSNE_PERPLEXITY=5

if [ ! -f ./FusionCatcher.tar.gz ]; then
    echo "Downloading test data set"
    wget_res=$(wget -q "https://dl.dropboxusercontent.com/u/16386922/co-fuse/FusionCatcher.tar.gz")
    if [ $? -ne 0 ]
    then
        echo "Fail to download the data set. The program will now terminate."
        exit 1;
    fi 
    tar xzf FusionCatcher.tar.gz
fi


if [ ! -d ./FusionCatcher ]; then
    echo "Folder 'FusionCatcher' cannot be found. Please unzip FusionCatcher.tar.gz"
    exit 1
fi


#
# Merge AML and MM results
#
mkdir -p FusionCatcher/AML_MM
cp -r FusionCatcher/AML/* FusionCatcher/AML_MM
cp -r FusionCatcher/MM/* FusionCatcher/AML_MM

#
# Run the script to reproduce the results in the paper
#

mkdir -p output_reproduce_results
Rscript --vanilla co-fuse.R 'FusionCatcher' './FusionCatcher/AML_MM' './output_reproduce_results/AML_MM' $TSNE_PERPLEXITY

Rscript --vanilla co-fuse2.R 'FusionCatcher' './FusionCatcher/AML' './FusionCatcher/MM' './output_fusioncatcher/fisher_test'

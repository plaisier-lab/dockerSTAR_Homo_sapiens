# dockerSTAR_Homo_sapiens
Dockerized STAR alignment tool that also builds the index for Homo sapiens GRCh38.p12.

# Ubuntu base
The latest clean ubuntu base is used as the base for the docker image

# Latest STAR installation
The latest STAR is downloaded from git as source code and compiled using the command:

```git clone https://github.com/alexdobin/STAR.git```

Source is compiled to executables and the location of the executables is added to the PATH.

# Homo sapiens GRCh38.p12 genome indexing
The files for the GRCh38.p12 genome are downloaded from [GENCODE](gencodegenes.org/human), ungzipped, and indexed into the folder `/index`.

# Build image using

```docker build -t star_2_7_1a_grch38_p21 .```

# Using the image with a mounted volume with your fastq files
```docker run -it -v /path/to/your/fastq/files:/fastq star_2_7_1a_grch38_p21 /bin/bash```

# Running STAR on the image
###TODO

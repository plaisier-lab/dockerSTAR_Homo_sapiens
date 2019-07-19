# dockerSTAR_Homo_sapiens
Dockerized STAR alignment tool that also builds the index for Homo sapiens GRCh38.p12.

# Ubuntu base
The latest clean ubuntu base is used as the base for the docker image

# Latest STAR installation
The latest STAR is downloaded from git as source code and compiled using the command:

'''git clone https://github.com/alexdobin/STAR.git'''

Source is compiled to executables and the location of the executables is added to the PATH.

# Homo sapiens GRCh38.p12 genome indexing
The files for the GRCh38.p12 genome are downloaded from [GENCODE](gencodegenes.org/human), ungzipped, and indexed into the folder `/index`.

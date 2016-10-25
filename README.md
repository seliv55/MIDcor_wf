![Logo](text4217.png)

# MIDcor
Version: 1.0
## “R”-program that corrects raw isotopic isomers spectra for natural occurring isotopes and overlapping of peaks for several metabolites in m/z scale

## Description

“midcor.R” is an “R”-program that performs a primary analysis of isotopic isomers (isotopomers) distribution obtained by Gas Cromatography coupled with Mass Spectrometry (GCMS). The aim of this analysis is to have a correct distribution of isotopes originated from substrates that are artificially enriched with specific isotopes (usually 13C). To this end the program performs a correction for natural occurring isotopes and also correction for “impurities” of the assay media that give peaks overlapping with the spectra of analyzed labeled metabolites. This program offers two ways of corrections of “impurities” resulted from overlapping the assayed mass isotopomer distribution with peaks produced either by unknown metabolites in the media, or by different fragments produced by the assayed metabolite. 

## Key features

- primary processing of 13C mass isotopomer data obtained with GCMS

## Functionality

- Preprocessing
- Statistical Analysis
- Workflows

## Approaches

- Isotopic Labelling Analysis
    - 13C
    
## Instrument Data Types

- MS


## Data Analysis

- correction for H+ loss produced by electron impact, natural occurring isotopes, and peaks overlapping

## Screenshots

- screenshot of input data (format Metabolights), output is the same format with one more column added: corrected mass spectrum

![screenshot](Screenshot.png)

## Tool Authors

- Vitaly Selivanov (Universitat de Barcelona)

## Container Contributors

- [Pablo Moreno](EBI)

## Website

- N/A

## Git Repository

- https://github.com/seliv55/midcor

## Installation

 # 1) As independent program. MIDcor itself does not require installation. Standing in the MIDcor directory enter in R environment with the command:
  
'''  R '''
  
  # read the necessary functions:
  
'''
source("lib.R")
  
source("midcor.R")
'''
  
  
     # 2) Docker image. To create the Docker container: i) go to the directory where the dockerfile is;
              ii) create container from dockerfile:
''' 
sudo docker build -t midcor:0.1 .
'''

## Usage Instructions

 # To run MIDcor independently: standing in the MIDcor directory inside R environment, after reading the sources execute the command:
 
 ''' run_midcor("input_file","output_file")  '''
 
 # here input file should be in Metabolights format, as is shown in the screenshot
 
 # To run MIDcor as a docker image, execute
 
 '''  sudo docker run -i -t -v $PWD:/data midcor:0.1 -i /data/input.csv -o /data/output.csv '''

 # An example of input file is provided as "outin.csv"

## Publications
- “MIDcor”, an R-program for deciphering mass interferences in mass spectra of metabolites enriched in stable isotopes. Submitted to BMC bioinformatics.

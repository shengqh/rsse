RNASeq Sample Size Estimation
==

* [Introduction](#Introduction)
* [Prerequisites](#Prerequisites)
* [Installation](#Installation)
* [Citation](#Citation)

<a name="Introduction"/>
#Introduction

Sample size calculation is an important issue in the experimental design
of biomedical research. For RNA-seq experiments, we already proposed a
method to calculate the sample size [[1](#Citation1)]. Here, we introduce our speed up
version of sample size estimation with user-friendly interface.

<a name="Prerequisites"/>
#Prerequisites

* Native library

GTK+ library is required to display the graphic interface. 

For linux system, the GTK+ usually has been already installed.

For windows system, user can download and intall GTK+ runtime enviroment from [https://sourceforge.net/projects/gtk-win/](https://sourceforge.net/projects/gtk-win/).

For mac system, it may be more diffcult to install correct GTK2 library. You may have a try using following command:
```
sudo port install gtk2 +x11
```

* R libraries
    * RGtk2, cairoDevice, gWidgets, and gWidgetsRGtk2 for user-friendly interface
    * stringr for parameter parsing
    * ssanv for RNAseq sample size estimation

You can install them by following command.

```
install.packages("gWidgetsRGtk2", dep = TRUE)
install.packages("stringr")
install.packages("ssanv")
```

<a name="Installation"/>
#Installation

RSSE package can be downloaded from [github](https://github.com/shengqh/rsse/releases)
Then, you can install and load RSSE package using following command:
```
install.packages(rsse_0.98.4.zip, repos = NULL)
library(rsse)
```

When you load the RSSE package on Windows system, it may ask you to install GTK+ libarary if there is no GTK libaray installed in your system. Please select GTK+ and follow the instruction to finish the GTK+ installation.

For mac system, it may be more diffcult to install correct GTK2 library.You may have a try on following command:
```
sudo port install gtk2 +x11
```

<a name="Citation"/>
#Citation
<a name="Citation1"/>
- Li CI, Su PF, Shyr Y: Sample size calculation based on exact test for assessing dierential expression analysis in RNA-seq data. BMC Bioinformatics 2013, 14:357.


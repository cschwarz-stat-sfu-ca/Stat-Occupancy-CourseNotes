Occupancy Estimation and Modeling.

The presence or absence of a species in a set of units (polygons, landscape units, territories,
etc.) is a fundamental concept in many ecological studies (e.g., resource selection modelling,
biodiversity, range). Visits to sampled units can result in a positive detection of a species
or a non-detection of the species in that unit. However, a species may not always be detected
if present which results in "false negatives". If the issue of detectablity is not accounted
for, estimates that rely on the observed (naive) level of occupancy can be misleading.

This course will cover methods for modelling species occurrence while accounting for potential
false negatives.


------------------------------------------------------------------------------

Please refer to
   http://www.cmiae.org/Events/#occupancy
for details on this course.

------------------------------------------------------------------------------

1. Introduction
   - differences/similarities of occupancy and capture-recapture modelling
   - fole of detectability
   - basic statistical review
   - concepts and notations
   - binomial data, odds ratios
   - basic likelihood
   - model selection and multi-model inference

2. Single-season occupancy models
   - basic sampling protocol
   - PRESENCE/MARK software to fitting single-season models
   - model assumptions
   - dealing with missing data
   - dealing with heterogeneity
   - more advanced features of PRESENCE/MARK - design matrices/covariates
   - planning studies; GENPRES software

3. Multiple-season occupancy models
   - sampling protocol
   - dynamics of population across seasons
   - using PRESENCE/MARK
   - alternative parameterizations
   - characterizing occupancy dynamics
   - modelling spatial correlations in occupancy dynamics
   - dealing with missing data
   - dealing with covariates
   - study design

4. Species co-occurrence models
   - sharing information among species
   - species richness or biodiversity
   - single season models/ multi season models
                                                                                                                                               
5. Extensions
   - modelling spatial correlations in occupancy dynamics
   - incorporation of count data
   - incorporation of marked animals
   - incorporation of telemetry data

------------------------------------------------------------------------------

Please download and install the following files. You may need admnistrator access
  to install the software.

  (a) Several software packages will be used.
      The following software is available for Windows operating systems only.

        Presence     https://www.mbr-pwrc.usgs.gov/software/presence.html
        GenPresence: https://www.mbr-pwrc.usgs.gov/software/presence.html
                     When you install Presence, GenPresence is automatically installed in the same director
                     as Presence

        MARK -       http://www.phidot.org/software/mark/downloads/

      The following software is available for Windoze and Mac operating systems

        RPresence    https://www.mbr-pwrc.usgs.gov/software/presence.html
        This is an R package that must then be installed in R.

        umarked      R package available from CRAN

        JAGS         http://www.stat.sfu.ca/~cschwarz/CourseNotes/HowGetSoftware.html
      Follow instructions at this website.


  (b) R/Rstudio      http://www.stat.sfu.ca/~cschwarz/CourseNotes/HowGetSoftware.html
      Also install the R packages as listed at this website.
      Be sure to create your own personal library as noted at the website.

      Install the RMark package from CRAN

      If you are using a Macintosh, it is possible to run RMark on a Mac directly
      without an emulator:

      # Refer to Simon Bonner instructions at:

         http://www.phidot.org/forum/viewtopic.php?f=21&t=3233#

      I had to change my shell to zsh before running Simon's scripts

     (1) Install home-brew
     /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

     (2) Install gcc
     brew install gcc

     (3) Install Mark
     brew tap sjbonner/tap
     brew install mark-on-mac

    (4) Check MARK - should return /usr/local/bin
    which mark

    (5) Check that mark runs - should return "no input file...."
    mark


    Everything seem to work now using RMark under Rstudio without any problems!
    Hurray,,, no more need for Windoze to run RMark! YMMV.



  (b) Adobe Reader: Course reading will be available in PDF format.
      You will need to call up some of the material to read from the
      screen during the course.
      You will be able to annotate the pdf files.

  (c) Excel: You will need to unzip files containing example datasets.
      Datasets are in an Excel workbook that you will need to open
      and “paste” into your program.

  (d) Text editor: Notepad, WordPad, or Word, to prepare data files for MARK/PRESENCE/GENPRES.

  (e) CourseMaterial  http://people.stat.sfu.ca/~cschwarz/CMIAE/
      You should download and unzip the file to get the full course material.

      You can also download individual files from the webpage.
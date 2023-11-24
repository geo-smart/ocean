



# 'Ocean': Understanding synoptic ocean sensor data


## Project Overview


The primary goal is to transform observational ocean data into
scientific insight into how the ocean works as a major component of the earth system.
Turning first to the data resource: The **Oceans Observatories Initiative** has operated 
since 2014 generating a profoundly rich dataset.
The variety and scope extend for nearly a decade across continental shelves and the deep ocean, across 
the water column from sea floor to surface, from seismology to marine mammal vocalizations to physical,
chemical and biological sensors. At the same time this complex system presents challenges in analysis 
and interpretation. 


The starting point or core of this project is a characterization of the ocean photic zone, from 200 meters
depth to the surface, the region of light. Regular observations include temperature, dissolved oxygen, 
salinity, particulate backscatter, chlorophyll concentration, dissolved organic matter fluorescence, 
nitrate concentration, pH, carbon dioxide concentration, spectral irradiance (7 wavelengths), 
spectrophotometry (2 forms x 83 channels), photosynthetically-available radiation and ocean current in 
three dimensions. These parameters are measured up to 18 times per day from 200 meters depth up to
near-surface in centimeter-scale increments. Consequently we are, roughly speaking, trying to 
understand the upper ocean in both depth and time as an 18-dimensional space. 


The interpretability of this 'starting point' synoptic dataset is itself daunting; but having come this
far we are also interested in making the problem "worse" after a fashion:


- Does the near surface time series agree with satellite observations?
- Does it agree with predictive ocean models?
- Does it agree with independent observations, for example from gliders or drifters (such as ARGO)?


Questions applicable to all of the data include


- How to get data in a usable form?
- Can the data be trusted?
- How can machine learning methods be constructively applied to this multi-dimensional space? 


## Elaboration: Potential Project Objectives


This project begins from the rich time series data described above (again: emphasizing the upper water column / photic zone). 
The following examples demonstrate potential project directions. Alternative ideas are welcome of course! 


1. Extend water column observations to map plane data using complementary observations from MODIS
2. Apply machine learning methods to *N > 3* coincident instrument data streams
3. Evaluate agreement between regional ocean model forecasts and *in situ* observation
4. Characterize water column biochemical profiles in relation to physical drivers such as turbulence and temperature
5. From real-time automated sensor data streams: Build a strategy for discrete sampling during maintenance cruises
6. Build an analysis workflow for automated sensors, e.g. anomaly/thin-layer detection
7. Use the data to differentiate influences of coastal upwilling versus terrigenous freshwater runoff


## Project Background: The Data Sketch


- Shallow profilers run a dozen or so instruments through the upper 200 meters of the water column nine times per day
- On the support platform at 200m depth an ADCP profiles current (200m to surface)
- ARGO floats: CTD and BGC
- ROMS ocean model
- MODIS sea surface data


## Project operational details


-


### Collaborators

-


### Specific questions / project goals


- specific tasks you want to accomplish, project goals, or research questions you want to answer.
- Think about what outcomes or deliverables you'd like to create
    - e.g. a series of tutorial notebooks demonstrating a [use case](https://geo-smart.github.io/usecases#Contributing)
    - e.g. a new python package


### Data


Briefly describe the data that will be used here (size, format, how to access).


### Existing methods


Needed: How would you or others traditionally try to address this problem?


### Proposed methods/tools


Needed: What new approaches would you like to try to implement to address your specific question(s) or application(s)?


### Additional resources or background reading


### Task Breakdown


### Repo structure

- Overview
    - Files
        - `.gitignore`
<br> Globally ignored files by `git` for the project.
        -  `environment.yml`
<br> `conda` environment description needed to run this project.
        - `README.md`
<br> Description of the project (see suggested headings below)

    - Folders
        - `contributors`
<br> Each team member can create their own folder under contributors, within which they can work on their own scripts, notebooks, and other files. Having a dedicated folder for each person helps to prevent conflicts when merging with the main branch. This is a good place for team members to start off exploring data and methods for the project.
        - `notebooks`
        - `scripts`
<br> Helper utilities that are shared with the team should go in here.
        - `data`
<br> Low-volume datasets that serve as excerpts from the full OOI operational time range. The sub-folders help
organize data by sensor/instrument/context/etcetera. 
        - `images`
<br> Resource images, for example invoked by notebook markdown for inline inclusion in a narrative.
        - `figures`
<br> Distinct from `images`: These are typically charts produced by notebook code execution.


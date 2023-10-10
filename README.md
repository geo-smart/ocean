



# 'Ocean': Understanding synoptic ocean sensor data


## Project Overview


The **Oceans Observatories Initiative** has operated since 2014 generating a profoundly rich dataset; and at the same
time a challenge in analysis and interpretation. The starting point is to characterize a single photic zone data set:
Regular observations of temperature, dissolved oxygen, salinity, particulate backscatter, chlorophyll 
concentration, dissolved organic matter fluorescence, nitrate concentration, pH, carbon dioxide concentration, 
spectral irradiance (7 wavelengths), spectrophotometry (2 x 83 channels), photosynthetically-available radiation
and ocean current in three dimensions: Measured 18 times per day over a depth range of the upper 200 meters of 
the ocean. This is roughly speaking depth/time-series in 18 dimensions. The interpretability of such a
synoptic dataset... rewrite paused here...


Questions as a means to making scientific progress:


- How are we to get data in a usable form?
- Can we trust the data once we have it?
- Once we have usable trustworthy data: What do we do with it? 


## *'The Problem'*: Potential Project Objectives


This project begins from a very rich set of data (emphasis: Upper water column / photic zone; see brief below) and provides five potential project directions to consider. Alternative ideas are welcome of course! The primary goal is to transform observational ocean data into
scientific insight into how the ocean works as a major component of the earth system.  


1. Extend water column observations to map plane data using complementary observations from the MODIS satellite
2. In search of synthesis: Apply machine learning methods to ten coincident instrument data streams
3. Evaluate agreement between regional ocean model forecasts and *in situ* observations
4. Characterize water column biochemical profiles in relation to physical drivers such as turbulence, salinity and temperature
5. Develop a discrete sampling strategy driven by real-time analysis of automated sensor data streams
6. Sharpen the analysis workflow for existing automated sensors
7. Investigate: Can the data be used to differentiate influences of coastal upwilling versus terrigenous freshwater runoff?

## Project Background: The Data Sketch


- und so weiter


## Project operational details


Should this project go forward: Here is where goes the details.


### Collaborators

List all participants on the project.

* Lead
* Team member
* Team member
* ...


### Specific questions / project goals


List the specific tasks you want to accomplish, project goals, or research questions you want to answer. Think about what outcomes or deliverables you'd like to create (e.g. a series of tutorial notebooks demonstrating a [use case](https://geo-smart.github.io/usecases#Contributing), or a new python package).


### Data


Briefly describe the data that will be used here (size, format, how to access).


### Existing methods


Needed: How would you or others traditionally try to address this problem?


### Proposed methods/tools


Needed: What new approaches would you like to try to implement to address your specific question(s) or application(s)?


### Additional resources or background reading

- und so weiter

### Task Breakdown


What are the individual tasks or steps that need to be taken to achieve the project goals? Think about which tasks are dependent on prior tasks, or which tasks can be performed in parallel. This can help divide up project work among team members.



- Repo structure
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


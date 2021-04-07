# SatellitePrediction [Version 1]
## Python 3.9
## Author: Nathan Aday / nraday1221@gmail.com
https://github.com/nathanaday/SatellitePrediction/

### DESCRIPTION
Generates upcoming satellite viewing opportunities at the user's location, with output observations in 
both a sqlite3 database and a .csv file. This program was created using Python 3.9

### USE
**To use SatellitePrediction, run Main.py and follow the prompts.**

For accurate results, you must get your own API key to access the elevation-api. It's FREE and takes only a minute.

- Step 1: Go to https://elevation-api.io/ and choose Signup/Login at the top

- Step 2: Create account with email and password (requires email verification; sometimes it says there is a problem 
confirming email, but you should be able to refresh and login without issue)

- Step 3: Create New API Key

- Step 4: Add new API key to **config.py** where it says:  *elevation_api = 'INSERT_API'*

- Step 5: Done, and now the program can access accurate elevation data for your location

*(Skipping this process will have the observation site altitude default to 0, which would impair accurate results 
at higher-altitude locations)*


**Make sure to satisfy package requirements in requirements.txt**

`pip install -r requirements.txt`





### DETAILED USE

There should be a directory in the project files called */output*. This is where the TLE satellite data, database, and 
obsfiles will be stored. The directory does not need to contain anything for the program to run, so it can be emptied 
as often as the user wants.

Further explanation of user inputs:

`Enter observation location (city, street, or nearby landmark...):`

The user can essentially enter anything that would get a hit on google, since it's using geopy with Nominatim 
open-street-map. The matched geocode location will appear (address, lat, lon, elv, and timezone). If it doesn't 
look right, the program should be run again, and the user should represent their location in a different context

`(Optional: Enter to skip) Name this site for the output file	    Ex: Home, Jakarta, Pikes Peak`
    
Whatever the user enters here is appended to the output .csv file. For the example of using "Home" as the site 
name, the user will see the file "OBS 04-06-2021(Home).csv" appear in /outputs when the program is complete.
This is especially helpful when creating observation files for multiple sites, but ultimately has no effect on the 
program's core functionality.


The easiest way to reference observations is to look at the .csv file. 


![Screen Shot 2021-04-06 at 6 09 10 PM](https://user-images.githubusercontent.com/79942554/113796069-3a13bc00-9703-11eb-9b3b-24c89793f12a.png)



Examplantion for output columns:

**Satellite Name,	Time,	Range(km),	Azimuth(deg),	Elevation(deg)**

ISS (ZARYA)	04/06/2021, 20:20:13,	1458.92,	299.76,	10.87

ISS (ZARYA)	04/06/2021, 20:21:53,	860.13,	283.32,	26.39

ISS (ZARYA)	04/06/2021, 20:23:33,	602.25,	214.77,	42.9

ISS (ZARYA)	04/06/2021, 20:25:13,	1000.67,	164.77,	20.96

- Column(1): The satellite name exactly as it appears in the TLE file.

- Column(2): The user’s local time at the moment the satellite has the shown range, azimuth, elevation values. Here we can 
expect the ISS to be visible for 5 minutes, between 20:20 and 20:25.

- Column(3): The range in kilometers between the observation site and the satellite.

- Column(4): The azimuth of the satellite’s location in the sky. 0 degrees is defined as north, so we can see here the ISS
will be transiting the sky from the Northwest (299) through true west (270) and toward the Southwest (164).

- Column(5): The elevation of the satellite’s location in the sky. The top of the user's visual hemisphere is defined as
90 degrees. In this example, we can see the elevation begins around 11 degrees and rises to a peak of about 43 degrees.


### PROGRAM CONTENTS
README.md 

requirements.txt 

config.py 

setup.py 

Data.py 

Main.py 

MathCore.py 

Satellite.py 

Site.py


### ADDITIONAL DETAILS
The satellite information comes from a two-line element (TLE) file compiled from NORAD data and hosted
on https://www.celestrak.com/NORAD/elements/. The TLEs on celestrack are updated several times a day, and the program
downloads the file at the beginning of every run. Therefore, the satellite data used within the program is always up to
date. For this particular program, the “100 (or so) Brightest” TLE is used, which includes anywhere between 80-200
satellites with a brightness (generally) sufficient for naked-eye observation.

The program propagates each satellite using 2-body EOM, considering both drag and J2 effects for orbital perturbation. 
A satellite is considered visible when it is not in the shadow of the Earth and it is above the user's horizon. Of 
course, magnitude of brightness plays a large role in naked-eye visibility, which is why the default TLE set comes from
the assorted bright list. There are still a variety of reasons why a satellite pass might not have good visibility, 
from weather to light-pollution.


### RESOURCES AND EXTENDED USE
The propagation methods used in this software are an application of the technqiues outlined in David A. Vallados’ 
Fundamentals of Astrodynamics and Applications, Four Edition.

For additional information on TLE and satellite tracking, as well as TLE files for a wide variety of satellite bodies,
see celestrack.com.

Changing what satellites this program propagates is as easy as changing the url in the download_tle() method located in
Data.py. Here are some URLs for quick reference:

Last 30 days launches: http://celestrak.com/NORAD/elements/tle-new.txt

Space Stations: http://celestrak.com/NORAD/elements/stations.txt

Starlink: http://celestrak.com/NORAD/elements/starlink.txt

### SUPPORT
Questions and bugs can be posted on the project's [github page](https://github.com/nathanaday/SatellitePrediction) or 
emailed to nraday1221@gmail.com

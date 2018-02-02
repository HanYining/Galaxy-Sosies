# Yining Han:
1. Sample link I know for now that gives a thumbnail of pictures. 
   http://skyserver.sdss.org/dr14/en/tools/chart/list.aspx
   
2. This should be the tool we use to do massive online batch query to our own
   Database. http://skyserver.sdss.org/CasJobs/
   
3. The guide on how to access the data base in batch mode:
   http://skyserver.sdss.org/CasJobs/Guide.aspx
   
   Check out this:
   The data base do have a nice documentation. (you will need to register)
   http://www.voservices.net/skyquery/Apps/Schema/
   
   I believe this is the main table we should use if we want to go with the images
   Table:  SDSSDR13:dbo.PhotoObjAll ( though i did not find a SDSSDR14 version)
   Primary Key: Unique SDSS identifier Flag: type = 3 to select galaxies.
   
4. Things I am trying to find out:
   a. what is the corresponding unique id for each galaxy.
   b. there is a query database limit 500MB, this should be a later concern.
   
5. I found this from the website:

   The SDSS project uses the color images to detect millions of objects (over 300 million stars and galaxies) over the whole northern sky. As soon as these objects are observed, a sophisticated software algorithm selects targets for further studies. Using another instrument, we take detailed spectra, measurements of the energy given off by the object as a function of its wavelength. An object's spectrum tells much more about an object than its image. A spectrum can be used to estimate a galaxy's distance and chemical composition, or to determine a star's age.
   
   So not sure if it is a good idea to start with image in the first place...


# Eva Gjekmarkaj

# Han Bao 

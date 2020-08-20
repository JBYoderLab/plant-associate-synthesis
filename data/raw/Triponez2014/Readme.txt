AFLP datasets
#############
STRUCTURE format, 
lines = specimens,
column 1 = specimen ID (three first capitals = population name)
column 2 = population number (for STRUCTURE use only)
columns >2 = AFLP data

files:
- Lysi439ind497loc.txt = L. vulgaris
- ME175ind51loc.txt = M. europaea
- MF74ind83loc.txt = M. fulvipes


Metadata
########
column 1 = population index
column 2 = population name (as coded in specimen IDs)
columns 3 - 5 = latitude, longitude and country

files:
- LysiMacro.info


Script
######
IBDpatterns_GST_GEO.r = pairwise GSTs among populations and relation with geographical distances


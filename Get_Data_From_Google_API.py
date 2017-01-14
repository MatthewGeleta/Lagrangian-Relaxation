#!python2
# Python 2 script for extracting Google Distance Matrix from list a town names

import urllib
import json
import scipy.io

# Google Distance Matrix API serivce URL:
serviceurl = 'https://maps.googleapis.com/maps/api/distancematrix/json?'
# Get list of locations from text file
fname = "Location_Names.txt"

fh = open(fname,'r')
# Put locations into Google Maps API format:
Num_Places = 0
Location_String = ""
for line in fh:
    Num_Places += 1
    Location_String += line.replace("\n", "|")
fh.close()

# Set string as list of origins and destinations
origins = Location_String
destinations = Location_String

# Set Google Distance Matrix parameters
mode = "walking" # eg: driving, bicycling, transit
language = "English"
units = "metric" # Note: only influences text field. Distances numerical values are stored as metric.
trafic_model = "best_guess" # eg: pessimistic, optimistic
transit_mode = "bus" # eg: subway, train, tram, rail or preference, train|tram|subway...
transit_route_preference = "fewer_transfers" # eg: less_walking
avoid = "" # eg: tolls, highways, ferries, indoor

# Encode information in url format
url = serviceurl + urllib.urlencode({'mode':mode,'sensor':'false', 'origins':origins, 'destinations':destinations})
print 
print 'Retrieving', url
connection = urllib.urlopen(url)
data = connection.read()
print 'Retrieved',len(data),'characters'

# Parse JSON to get python dictionary
try: js = json.loads(str(data))
except: js = None
if 'status' not in js or js['status'] != 'OK':
    print '==== Failure To Retrieve ===='
    print data
    exit()
    
# Number of retrieved locations:
Num_Places = len(js["rows"])
# Initialise distance matrix (square matrix, entries are travel times)
Dist_Matrix = [[0 for x in range(Num_Places)] for y in range(Num_Places)]

# Retrieve travel times (in seconds) between all retrieved locations:
for row in range(0, Num_Places):
    for col in range(0, Num_Places):
        duration = js["rows"][row]["elements"][col]["duration"]["value"]
        Dist_Matrix[row][col] = duration
        
# Save data to MATLAB format
scipy.io.savemat("MATLAB_Info", {"Dist_Matrix": Dist_Matrix, "Num_Places": Num_Places})
        
# Uncomment the next two lines to pretty-print all of the retrieved JSON
#print "Full GoogleMaps info:"
#print json.dumps(js, indent=4)

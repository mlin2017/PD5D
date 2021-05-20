import os
import re
for filename in os.listdir("."):
    if re.match("annotations_*", filename):
        thisfile=open(filename)
        testout=open("aggrtable.csv","w")
        next(thisfile)
        testout.write("library_id,molecule_h5\n")
        for line in thisfile:
            tmpline=(line.split(",")[1],'../cellranger_count/%s/outs/molecule_info.h5\n' % line.split(",")[1])
            fnlline=",".join(tmpline)
            testout.write(fnlline)

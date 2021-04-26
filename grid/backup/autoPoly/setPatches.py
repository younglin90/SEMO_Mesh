### script to determine liquidus and solidus position over time

import os

boundaryfile = open("./constant/polyMesh/boundary", 'r')
content = boundaryfile.read()
boundaryfile.close()

#fix front
pos = content.find("front")
print "POS1: "+ str(pos)
pos2 = content.find("patch", pos, pos + 50)
print "POS2: "+ str(pos2)
#sanity check - file could have already been processed
if (pos2 == -1):
	print "Patch not found - quitting"
	quit()
content = content[0:pos2] + "empty" + content[pos2+5:]

#fix back
pos = content.find("back")
print "POS1: "+ str(pos)
pos2 = content.find("patch", pos, pos + 50)
print "POS2: "+ str(pos2)
content = content[0:pos2] + "empty" + content[pos2+5:]

#fix defaultFaces
pos = content.find("defaultFaces")
print "POS1: "+ str(pos)
pos2 = content.find("patch", pos, pos + 50)
print "POS2: "+ str(pos2)
content = content[0:pos2] + "empty" + content[pos2+5:]

boundaryfile = open("./constant/polyMesh/boundary", 'w')
boundaryfile.write(content)
boundaryfile.close()

print "DONE."

quit()



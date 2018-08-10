#!/usr/bin/python3

import dbm
import os
import sys

path=os.getcwd()
#files = os.listdir(path)
accFile = sys.argv[1]

with dbm.open('rev_accDB.db','c') as db:
#    for x in files:
#        if x.endswith(".txt"):
#            filePath = path + "/" + x
#            GCA = x[0:15]
#            with open(filePath,"r") as y:
#                for line in y:
#                    db[line] = GCA
#
#            db[GCA]=genomeVersion
#    print(db.get(b'GCA_000513895.1').decode('utf-8'))
	with open(accFile, "r") as y:
		for line in y:
			try:
				print(line.rstrip() + "\t",db.get(line).decode('utf-8'))
			except:
				pass
		

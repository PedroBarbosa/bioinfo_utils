#!/usr/bin/env python

import sys
from sys import version_info
from sys import stdin

import exceptions
import time
import argparse

import struct
from struct import *

try:
	libreadline = __import__("readline")
except ImportError as e:
	pass


try:
	liblmdb = __import__("lmdb")
except ImportError as e:
	print "Error: %s" % e
	print "failed to load lmdb module. Install it using \"$ pip install lmdb\" command\n"
	exit(1)
else:
	globals()["lmdb"] = liblmdb


parser = argparse.ArgumentParser(description='lookup for accession.', version="%s version %s" % (sys.argv[0], "1.1"))
parser.add_argument('-d', '--database', dest='dbpath', default='./gi2acc_lmdb.db', help='path to database file, ./gi2acc_lmdb.db is used by default')
parser.add_argument('-i', '--info', dest='printinfo', action="store_true", help='print meta data from db')
args = parser.parse_args()
py3 = version_info[0] > 2

filename=args.dbpath
GI_DBI='#GI_DBI'
META_DBI='#META_DBI'
META_INCREMENTAL_TIME='INC_TIME'
META_FULL_TIME='FULL_TIME'
if stdin.isatty():
	PROMPT='gi: '
else:
	PROMPT=''

try:
	env = lmdb.open(filename, subdir=False, readonly=True, create=False, max_dbs=16)
	gidb = env.open_db(GI_DBI, create=False)
	if args.printinfo:
		print "Opened %s" % filename
except Exception as e:
	print "Failed to open %s\nError: %s" % (filename, e)
	exit(1)

try:
	metadb = env.open_db(META_DBI, create=False)
	txn = env.begin(metadb)
	tmfull = txn.get(key = META_FULL_TIME, db=metadb)
	tminc = txn.get(key = META_INCREMENTAL_TIME, db=metadb)
	txn.commit();
	if args.printinfo:
		s = "Created: " + time.ctime(int(tmfull))
		if tminc:
			s += ", last updated: " + time.ctime(int(tminc))
		print "%s\n" % s
		env.close()
		exit(0)

except Exception as e:
	pass


while True:
	try:
		if py3:
			s = input(PROMPT)
		else:
			s = raw_input(PROMPT)
		if not s : continue
		if s[0] == 'q' : break
		gi = int(s)
	except KeyboardInterrupt:
		print "\nfinished\n"
		exit(0)
	except EOFError:
		print "\n"
		exit(0);
	except Exception as e:
		print "Error: %s" % e
		exit(1)

	try:
		gikey = pack('<q', int(gi))
#		print ":".join('{:02x}'.format(ord(c)) for c in gikey)
		txn = env.begin(gidb)
		buf = txn.get(key=gikey, db=gidb)
		txn.commit();
		if buf:
#			print ":".join('{:02x}'.format(ord(c)) for c in buf)
			gi_len_flg = ord(buf[0])
			gi_len_bytes = gi_len_flg & 0x7
			gi_len_isneg = gi_len_flg & 0x8
			gi_len = 0
			if (gi_len_bytes > 0):
				shift = 0
				for n in buf[1:1+gi_len_bytes]:
					gi_len += ord(n) << shift
					shift += 8
			if (gi_len == 0) and gi_len_isneg:
				gi_len = -1
			elif gi_len_isneg:
				gi_len = -gi_len;
				
			acc_len = ord(buf[1 + gi_len_bytes])
			acc_ofs = 2 + gi_len_bytes
			acc = buf[acc_ofs:acc_ofs + acc_len]
			print "%d\t%s\t%d" % (gi, acc, gi_len)
		else:
			print "%d not found" % (gi)

	except Exception as e:
		print "Error: %s" % e


env.close()



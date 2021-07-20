# used by HLAminer.py
import sys
import os

cwdir=sys.argv[1]
cmdline=sys.argv[2]

os.chdir(cwdir)
os.system('sh ' + cmdline)

#!/usr/bin/python

'''
Agent configuration tool; this checks cc is available and adds the fileshare to the sandbox file,
using the configuration file available in the share
'''

import HAConf
import os

shared_location = HAConf.sidekick_location[:-8]

print "Configuring xgrid sandbox"
print "I'm going to ask your password now to change a system file relating to xgrid (using sudo)"

if shared_location[-1] != "/":
	additional_xgrid_line = '(allow file-read* file-write* (regex "^'+shared_location+'(/|$)")) ; Sidekick'
else:
	additional_xgrid_line = '(allow file-read* file-write* (regex "^'+shared_location[:-1]+'(/|$)")) ; Sidekick'
	
os.system("""(grep -v Sidekick /usr/share/sandbox/xgridagentd_task_nobody.sb ; echo '""" + additional_xgrid_line + """') > /tmp/xgridagentd_task_nobody.sb; sudo mv /tmp/xgridagentd_task_nobody.sb /usr/share/sandbox/xgridagentd_task_nobody.sb""")


#Finally, run some tests to see if we have everything we need

try:
	import matplotlib, pylab
except:
	print "WARNING:	Matplotlib is not installed; graphics cannot be generated"

if not os.path.exists("/usr/bin/cc"):
	print "WARNING:	Xcode has not been installed; gromacs can't run without /usr/bin/cc" 
	
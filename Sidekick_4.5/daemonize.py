import sys,os

def daemonize (stdin='/dev/null',stdout='/dev/null', stderr='/dev/null'):
	#perform first fork
	try:
		pid = os.fork()
		if pid > 0:
			sys.exit(0) #Exit first parent
	except OSError, e:
		sys.stderr.write("Fork 1 failed: (%d) %s\n" % (e.errno, e.strerror))
		sys.exit(1)
	#Decouple from parent environment
	os.chdir("/")
	os.umask(0)
	os.setsid()
	#Perform second fork
	try:
		pid = os.fork()
		if pid > 0:
			sys.exit(0) #Exit first parent
	except OSError, e:
		sys.stderr.write("Fork 2 failed: (%d) %s\n" % (e.errno, e.strerror))
		sys.exit(1)
	#The process is now daemonized, redirect standard file descriptors
	for f in sys.stdout, sys.stderr: f.flush()
	si = file(stdin, 'r')
	so = file(stdout, 'a+')
	se = file(stderr, 'a+',0)
	os.dup2(si.fileno(),sys.stdin.fileno())
	os.dup2(so.fileno(),sys.stdout.fileno())
	os.dup2(se.fileno(),sys.stderr.fileno())

def _example_main():
	import time
	sys.stdout.write("Daemon started with pid %d\n" % os.getpid())
	sys.stdout.write("Daemon stdout output\n")
	sys.stderr.write("Daemon stderr output\n")
	c = 0
	while True:
		sys.stdout.write('%d: %s\n' % (c, time.ctime()))
		sys.stdout.flush()
		c = c+1
		time.sleep(1)

if __name__ == "__main__":
	daemonize('/dev/null','/tmp/daemon.log','/tmp/daemon.log')
	_example_main()
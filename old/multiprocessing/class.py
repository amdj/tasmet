import multiprocessing
import time
class Worker(multiprocessing.Process):
	def __init__(s,i):
		multiprocessing.Process.__init__(s)
		s.i=i
	def run(s):
		time.sleep(s.i)
		print 'In %s' % s.name
		return

if __name__ == '__main__':
	jobs = []
	for i in range(5):
		p = Worker(5-i)
		jobs.append(p)
		p.start()
	for j in jobs:
		j.join()

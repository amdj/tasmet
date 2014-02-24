import multiprocessing as mp

class ToBeUsed:
	def __init__(s,name):
		s.name=name
	def do_something(s,dataneeded):
		proc_name=mp.current_process().name
		print "Process name: %s" %proc_name
		print "Object name: %s" %s.name
		print "Data: %g" %data

def worker(queue,data):
	objec=queue.get()
	objec.do_something(data)

if __name__=="__main__":
	queue=mp.Queue()
	data=1
	p=mp.Process(target=worker, args=(queue,data))
	p.start()
	queue.put(ToBeUsed('test'))
	queue.close()
	queue.join_thread()
	p.join()

import multiprocessing
import time

class Consumer(multiprocessing.Process):
	def __init__(self, task_queue, result_queue):
		multiprocessing.Process.__init__(self)
		self.task_queue = task_queue
		self.result_queue = result_queue

	def run(self):
		proc_name = self.name
		while True:
			next_task = self.task_queue.get()
			if next_task is None:
			# Poison pill means we should exit
				print '%s: Exiting' % proc_name
				break
			print '%s: %s' % (proc_name, next_task)
			answer = next_task()
			self.result_queue.put(answer)
		return

class Task:
	def __init__(self, a, b):
		self.a = a
		self.b = b
	def __call__(self):
		time.sleep(0.1) # pretend to take some time to do our work
		return '%s * %s = %s' % (self.a, self.b, self.a * self.b)
	def __str__(self):
		return '%s * %s' % (self.a, self.b)


if __name__ == '__main__':
	# Establish communication queues
	tasks = multiprocessing.Queue()
	results = multiprocessing.Queue()
	# Start consumers
	num_consumers = 4#multiprocessing.cpu_count() * 2
	print 'Creating %d consumers' % num_consumers
	consumers = [Consumer(tasks, results) for i in xrange(num_consumers) ]
	#^Creates an array of 4 consumers, give each consumer the task queue
	# and the results queue
	for w in consumers:
		w.start()
		# Enqueue jobs
	num_jobs = 10
	for i in xrange(num_jobs):
		tasks.put(Task(i, i))
	# Add a poison pill for each consumer
	for i in xrange(num_consumers):
		tasks.put(None)

	# Start printing results
	while num_jobs:
		result = results.get()
		print 'Result:', result
		num_jobs -= 1


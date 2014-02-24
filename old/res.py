'''
Created on Jul 30, 2013

@author: anne
'''
import numpy as n

class system(object):
	def __init__(self,segments):
		self.nsegments=len(segments)
		self.segments=segments
		self.dofs=[]
		self.totdofs=0 #Total number of dofs
		self.cumdofs=[0] #Cumulative dofs
		j=0
		self._bc=[]
		number0=segments[0].number

		segnum=0
		for segment in self.segments:
			segment.number=segnum
			self.dofs.append(segment.ndofs)
			self.totdofs+=segment.ndofs
			self.cumdofs.append(self.cumdofs[j-1]+segment.ndofs)

		self._result=n.zeros((self.totdofs),float)
		self._error=n.ones((self.totdofs),float)
		for segment in self.segments:
			segnum=segment.number
			start,end=self.getSegRange(segnum)
			segres=self.getSegResult(segnum)
			self._result[start:end]=segres
			segnum+=1

	def deleteSegment(self,number):
		del self.segments[number]
		j=0
		for i in self.segments:
			i.number=j
			j+=1
		del self._bc
		self._bc=[]
	def reset(self):
		for i in self.segments:
			i.reset()
		self.__init__(self,self.segments)
	def getSegRange(self,segnum):
		start=self.cumdofs[segnum]
		end=start+self.dofs[segnum]
		return start,end
	def getSegResult(self,segnum):
		start,end=self.getSegRange(segnum)
		return self._result[start:end]
	def setSegResult(self,segnum,segres):
		start,end=self.getSegRange(segnum)
		self._result[start:end]=segres
		self.segments[segnum].setResult(segres)
	def setResult(self,res):
		#Here, we can paralellize
		self._result=res
		for segment in self.segments:
			segnum=segment.number
			start,end=self.getSegRange(segnum)
			segres=res[start:end]
			self.setSegResult(segnum, segres)
	def getResult(self):
		return s._result[:]
	def getSegError(self,segnum):
		self.segments[segnum].setResult(segresult)
		return self.segments[segnum].getError()
	def error(self,res):
		#This can, again be parallelized
		self.setResult(res)
		for i in xrange(self.nsegments):
			segerror=self.getSegError(i)
			start,end=self.getSegRange(i)
			self._error[start:end]=segerror
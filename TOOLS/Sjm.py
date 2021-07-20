# coding=utf-8
# tools for SJM jobs construction/parse/list/write
import sys
import os.path
import os

class job:   #SJM 任务 类
	def __init__(self, name, cmd, memory, status='waiting', sched_options='-V -cwd', cpu = 1, q='all.q', path = None):
		self.name=name
		self.memory=memory
		self.cpu = cpu
		self.status=status
		self.sched_options = sched_options+' -q ' + q + ' -l p=' + str(self.cpu) 
		self.cmd=cmd
		self._dataValidation()
		self._shellMake()

	def _shellMake(self):
		if path:
			if not os.path.exist(path):
				os.mkdir(path)
			jobshellFile = os.path.join(path, self.name + '.sh')
			with open(jobshellFile, 'w') as jobShell:
				jobShell.write(self.cmd)
			self.cmd = 'sh ' + jobshellFile

	def _dataValidation(self):
		if self.memory[-1].upper() not in ['G', 'M']:
			sys.exit('illegal memory setup ...')
		
		
	def printJob(self):
		code='''
job_begin
  name %s
  memory %s
  status %s
  sched_options %s
  cmd_begin
    %s
  cmd_end
job_end
''' % (self.name, self.memory, self.status, self.sched_options, self.cmd)
		return code


class order:   #SJM 顺序 类
	def __init__(self, jobA, jobB, type='before'):
		self.jobA=jobA
		self.jobB=jobB
		self.type=type
		self._dataValidation()
		
	def _dataValidation(self):
		if self.type not in ['before', 'after']:
			sys.exit('[ERROR] illegal order keyword ...')
			
	def showOrder(self):
		if type == 'before':
			return [self.jobA, self.jobB]
		else:
			return [self.jobB, self.jobA]
			
	def printOrder(self):
		return ' before '.join(self.showOrder())+'\n'
			
			
class sjmUnit:   #SJM任务单元（job & order） 类
	def __init__(self, job, order = []):
		self.job=job
		self.order=order
		self._dataValidation()
		
	def _dataValidation(self):
		if not isinstance(self.job, job):
			sys.exit('[ERROR] illegal object job for sjmUnit ...')
		for eachOrder in self.order:
			if not isinstance(eachOrder, order):
				sys.exit('[ERROR] illegal object order for sjmUnit')
			if self.job.name not in eachOrder.showOrder():
				sys.exit('[ERROR] Non-match job and order object in sjmUnit ...')
				

class sjmFile(list):   #SJM job文件 类  继承list类
	def _dataValidation(self):
		nameFromJob=[]
		nameFromOrder=[]
		orderRule=[]
		for eachUnit in self:
			assert isinstance(eachUnit, sjmUnit)   #确认是正确的类
			nameFromJob.append(eachUnit.job.name)
			for eachOrder in eachUnit.order:
				nameFromOrder+=eachOrder.showOrder()
				orderRule.append(eachOrder.showOrder())
		if not set(nameFromOrder) == set(nameFromJob):   #不能要求完全相等，还要衔接其他的流程一起的
			print '[WARNING] Non-match job and order object in SJM file ...'
			
		nameDedup=[]
		for eachName in nameFromOrder:   #保持原定顺序，所以没有使用set  两次遍历法做判断 1. 构建顺序序列  2. 再次判读矛盾order
			if eachName not in nameDedup:
				nameDedup.append(eachName)	
		for eachRule in orderRule:
			index0=nameDedup.index(eachRule[0])
			index1=nameDedup.index(eachRule[1])
			if  index0 > index1:
				(nameDedup[index0], nameDedup[index1]) = (nameDedup[index1], nameDedup[index0])
		for eachRule in orderRule:
			index0=nameDedup.index(eachRule[0])
			index1=nameDedup.index(eachRule[1])
			if index0 > index1:
				sys.exit('[ERROR] conflicting job order ...')

	def printSjmJob(self):
		self._dataValidation()
		jobs=''
		orders=''
		for eachUnit in self:
			jobs+=eachUnit.job.printJob()
			for eachOrder in eachUnit.order:
				orders+=eachOrder.printOrder()
		return '\n'.join([jobs,orders])
		
		

def parseJobStatus(file):
	if not os.path.isfile(file):
		sys.exit('[ERROR] ' + file+' not exist ...')
	jobStatusDict={}
	logDir=''
	with open(file) as statusFile:
		flag=0 
		for eachLine in statusFile:
			if eachLine.startswith('log_dir'):
				logDir=eachLine.strip().split()[1]
				continue	
			if eachLine.startswith('job_begin'):
				flag=1 
				jobName=''
				jobStatus=''
			if flag and eachLine.strip().startswith('name'):
				jobName=eachLine.strip().split()[1]
			if flag and eachLine.strip().startswith('status'):
				jobStatus=eachLine.strip().split()[1]
			if eachLine.startswith('job_end'):
				flag=0 
				if jobName not in jobStatusDict:
					jobStatusDict[jobName]=jobStatus
	return logDir, jobStatusDict



if __name__ == '__main__':
	jobdemo=job('jobA', '500M', 'perl try.pl')
	jobdemo2=job('jobB', '500M', 'perl try2.pl')
	orderdemo=order('jobA', 'jobB', 'after')
	unitdemo1=sjmUnit(jobdemo)
	unitdemo=sjmUnit(jobdemo2,orderdemo)
	f=sjmFile()
	f.append(unitdemo)
	f.append(unitdemo1)
	print f.printSjmJob()
	
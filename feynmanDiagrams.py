from __future__ import division
from guppy import hpy
import particlePhysics
import rambo
import numpy as np
import random
import time
import gc
import operator
import SED


VERTEX_ID = 0
PRINT = False
PRINT_ALL = False

# Holds multiple lines
class Vertex():
	def __init__(self):
		self.neighbours = []
		global VERTEX_ID
		VERTEX_ID += 1
		self.vertexID = VERTEX_ID
	def __str__(self):
		return "Vertex%s" %self.vertexID
	def __repr__(self):
		return self.__str__()
	def addLine(self, line, neighbour=None):
		if len(self.neighbours) < 3:
			if neighbour == None:
				self.neighbours.append([neighbour, line])
				return
			elif len(neighbour.neighbours) < 3:
				self.neighbours.append([neighbour, line])
				neighbour.neighbours.append([self, line.getReversed()])
				return
		print "ERR: Too many lines at one vertex!!!"
				
				
# Holds particle and it's direction		
class Line():
	def __init__(self, particle, direction):
		self.particle = particle
		self.direction = direction
		if particle.onShell:
			self.cacheCode = int(2**particle.particleNum)
	def __str__(self):
		return "%s (dir:%s)" %(self.particle.name, self.direction)
	def __repr__(self):
		return self.__str__()
	def getReversed(self):
		return Line(self.particle, -1*self.direction)
	def getNormalizedCacheCode(self, maximum):
		if self.cacheCode >= maximum - self.cacheCode: return maximum - self.cacheCode
		else: return self.cacheCode
		
# Holds list of vertices and propagator binary codes	
class FeynmanDiagram():
	def __init__(self, vertices, propagators):
		self.vertices = vertices
		self.propagators = propagators
		self.matrixElement = 1
		self.size = 0
	# Traverse diagram, printing out structure (for debugging)
	def traverse(self):
		for vertex in self.vertices:
			print "\nVertex:"
			for line in vertex:
				particle = line.particle
				if line.direction == -1: print "%s is going out" %(particle.name)
				else: print "%s is coming in" %(particle.name)
		print "Matrix Element: %s" %(self.matrixElement)

	def produceList(self):
		vertexLists = []
		currentList = []
		for vertex in self.vertices:
			cacheCode = 0
			for line in vertex[:-1]:
				cacheCode += line.cacheCode 
			currentList = currentList + [cacheCode]
			vertexLists.append(tuple(currentList))
		return vertexLists
		
	def connectUp(self):
		for vertex in self.vertices:
			for line in vertex:
				particle = line.particle
				if line.direction == -1:
					particle.goingToVertex = vertex
				else:
					particle.comingFromVertex = vertex	

# Builds diagrams and calculates matrix elements
class Interaction():	
	
	def __init__(self, incomingParticles, outgoingParticles, model):
		self.lineList = []
		self.matrixElement = 0
		self.incomingEnergy = 0
		#model.onNewInteraction(incomingParticles, outgoingParticles)
		self.model = model
		#self.incomingMom = particlePhysics.ThreeVector(0, 0, 0)
		numP = 0
		for particle in incomingParticles:
			if particle.name[-1] != ")": particle.name += " (%s)" %(numP+1)
			particle.particleNum = numP
			self.lineList.append( Line(particle, 1) )
			self.incomingEnergy += particle.energy
			numP += 1
		#~ self.sijThreshold = self.incomingEnergy**2 - self.incomingMom.dot(self.incomingMom)
		for particle in outgoingParticles:
			if particle.name[-1] != ")": particle.name += " (%s)" %(numP+1)
			particle.particleNum = numP
			self.lineList.append( Line(particle, -1) )
			numP += 1
		self.cacheMax = 0
		for j in range(numP): self.cacheMax += 2**j
		self.cache = [None for i in range(self.cacheMax+1)]
		for line in self.lineList:
			self.addToCache(line)
		self.particleCache = [[] for i in range(self.cacheMax+1)]
		self.sortingTime = 0
		self.diagCalcTime = 0
		self.particleCalcTime = 0
		self.recursionTime = 0
		self.lastRecursionTime = 0
		
		self.diagramSet    = set()
		self.diagramSets   = set()
		self.diagramSetLen = 0
		self.diagramList   = []
		self.totalDiagsBuilt = 0
		self.matrixElement = 0
		
		self.faster = True
		
	def isAllowed(self):
		charge = 0
		energy = 0
		momentum = particlePhysics.ThreeVector(0,0,0)
		threshold = self.incomingEnergy * 1e-11
		for line in self.lineList:
			direction = line.direction
			particle  = line.particle
			#print particle.energy, particle.threeMomentum
			energy += particle.energy * direction
			charge += particle.charge * direction
			momentum.x += particle.threeMomentum.x * direction
			momentum.y += particle.threeMomentum.y * direction
			momentum.z += particle.threeMomentum.z * direction
		if abs(energy) > threshold:
			return False
		if momentum.mod() > threshold:
			return False
		if charge != 0:
			return False
		return True
		
	def changed(self):
		self.diagramSet    = set()
		self.diagramSetLen = 0
		self.diagramList   = []
		self.totalDiagsBuilt = 0
		self.matrixElement = 0
		self.diagramSets   = set()
	
	def calcDiagrams(self, saveVertices = False, optimize = False, diagNumber=False):
		self.diagNumber = diagNumber
		self.buildDiagrams(self.lineList, first=True, saveVertices = saveVertices, optimize=optimize) 
		self.diagrams = self.diagramList	
		return self.diagrams
		
	def getDiagrams(self):
		self.buildDiagrams(self.lineList, first=True, saveVertices = True) 
		self.diagrams = self.diagramList	
		return self.diagrams
		
	def isDuplicate(self, line):
		if self.cache[line.cacheCode] != None:
			if self.cache[line.cacheCode].direction != line.direction:
				return True
		return False
		
	def getFromCache(self, cacheCode):
		line = self.cache[cacheCode]
		return line
		
	def getNormalizedCacheCode(self, cacheCode):
		cache2 = self.cacheMax - cacheCode
		if cacheCode > cache2:
			return cache2
		else:
			return cacheCode
	
	def addToCache(self, line):
		otherCode = self.cacheMax - line.cacheCode
		line2 = self.cache[line.cacheCode]
		if line2 == None: # Need to add the reverse in the -ve dir.
			self.cache[line.cacheCode] = line
			self.cache[otherCode] = line.getReversed()
			self.cache[otherCode].direction = line.direction
			self.cache[otherCode].particle = particlePhysics.MasslessScalar(onShell=False)
			self.cache[otherCode].particle.energy = -line.particle.energy
			self.cache[otherCode].particle.threeMomentum = line.particle.threeMomentum.getScaled(-1)
			self.cache[otherCode].cacheCode = otherCode
		elif PRINT_ALL:
			print "Err, line already exists in cache"
			
	def addToParticleCache(self, line):
		otherCode = self.cacheMax - line.cacheCode
		line2 = self.cache[line.cacheCode]
		if line2 == None:
			self.cache[line.cacheCode].app
		

	def buildDiagrams(self, lineList, first=False, propagator = 1, propagators = [], rawPropagators=[], vertices=[], fourVertices = 0, saveVertices=False, optimize=False):
		if first:
			if not self.isAllowed():
				print "Err: Conservation Laws not obeyed, interaction impossible"
				return
		if PRINT: 
			"Building started"
		if first:
			startTime = time.time();
		# If only 3 or 4 lines, try and create a working vertex
		if len(lineList) == 3 or (len(lineList) == 4 and self.model.fourVertices):
			if (PRINT_ALL): print "success"
			debug = False
			newPropagator = self.model.vertexPropagator([ [line.particle, line.direction] for line in lineList ], None, propagator)
			if newPropagator != None and newPropagator != "Err": 
				if saveVertices: self.totalDiagsBuilt += 1
				propagators.sort()
				traverse = False
				if newPropagator.real == 0: newPropagator = 0 + newPropagator.imag*1j
				if newPropagator.imag == 0: newPropagator = newPropagator.real + 0*1j
				roundedPropagator = '%.4g + %.4g j' %(newPropagator.real, newPropagator.imag)
				self.diagramSet.add(tuple(propagators + [roundedPropagator]))
				length = len(self.diagramSet)
				if  length > self.diagramSetLen:
					self.diagramSetLen = length 
					if saveVertices: 
						diagram = FeynmanDiagram(vertices+[lineList], propagators)
						if len(lineList) == 4: 
							diagram.fourVertices = fourVertices + 1
						else:
							diagram.fourVertices = fourVertices
						diagram.matrixElement = newPropagator
						self.diagramList.append(diagram)
					self.matrixElement += newPropagator
					if optimize and self.diagNumber != False:
						if length == self.diagNumber:
							return False
			# Don't return because still need to do the below
		# If more than 3 lines
		if len(lineList) > 3:
			diagrams = []
			tempDiagrams = []
			# Cycle over possible line pairs
			for index1, line1 in enumerate(lineList): 
				if first and PRINT and not PRINT_ALL: print "%s of %s" %(index1, len(lineList))
				for index2, line2 in enumerate(lineList):
					if ( index2 <= index1 ): 
						continue
					if PRINT_ALL and first: print "Top-level"
					newCacheCode = line1.cacheCode + line2.cacheCode
					newCacheCodeNormed = self.getNormalizedCacheCode(newCacheCode)
					newPropList = propagators+[newCacheCodeNormed]
					newPropList.sort()
					if optimize == False:
						newRawPropList = None
					else:
						newRawPropList = rawPropagators + [newCacheCode]
					if optimize == False:
						nVert = len(newPropList)
						setLen = len(self.diagramSets)
						self.diagramSets.add(tuple(newPropList+[fourVertices]))
						threeVert = True
						if self.faster and len(self.diagramSets) == setLen:
							threeVert = False
					else:
						if tuple(newRawPropList) in optimize:
							threeVert = True
						else:
							threeVert = False
					particleStartTime = time.time()
					particle1 = line1.particle
					particle2 = line2.particle
					direction1 = line1.direction
					direction2 = line2.direction
					if (particle1.particleType in particle2.couplesWith) and (particle2.particleType in particle2.couplesWith): # Make sure the 2 particles can couple
						if threeVert:
							# Try and get the line from the cache
							placeholderParticle = particlePhysics.MasslessScalar(onShell=False)
							placeholderLine = Line( placeholderParticle, 1 )
							placeholderLine.cacheCode = newCacheCode
							cacheLine = self.cache[newCacheCode]
							# If it's not there, calculate it and add it
							if cacheLine == None:		
								placeholderParticle.threeMomentum.x = (particle1.threeMomentum.x*direction1 + particle2.threeMomentum.x*direction2) 
								placeholderParticle.threeMomentum.y = (particle1.threeMomentum.y*direction1 + particle2.threeMomentum.y*direction2) 
								placeholderParticle.threeMomentum.z = (particle1.threeMomentum.z*direction1 + particle2.threeMomentum.z*direction2) 
								placeholderParticle.energy = (particle1.energy * direction1 + particle2.energy * direction2)
								self.addToCache(placeholderLine)
							# If it is, copy it.
							else:
								placeholderParticle.threeMomentum = cacheLine.particle.getThreeMomentum()
								placeholderParticle.energy = cacheLine.particle.energy
							# Get the 3rd particle type
							listNewLines = []
							vertexCharge = particle1.charge*direction1+particle2.charge*direction2
							for particleType in particle1.couplesWith.intersection(particle2.couplesWith):
								# Can't have a vertex with all particles the same particle type (unless it's a Phi particle)
								if particleType != particlePhysics.Phi and ( (particle1.particleType == particleType) and (particle2.particleType == particleType) ):
									self.particleCalcTime += time.time() - particleStartTime
									continue
								# Creating new particle, and a new incoming line for the rest of the diagram
								newParticle = particleType(onShell=False) # This is where I need to be checking the momentum cache instead of calculating
								newParticle.threeMomentum = placeholderParticle.threeMomentum
								newParticle.energy = placeholderParticle.energy
								# Check conservation laws are obeyed
								if newParticle.charge != vertexCharge:
									self.particleCalcTime += time.time() - particleStartTime
									continue	
								if (PRINT_ALL): print newLineLists
								newLine = Line( newParticle, 1)
								newLine.cacheCode = newCacheCode
								listNewLines.append(newLine)
							for newLine in listNewLines:
								newParticle = newLine.particle
								newPropagator = self.model.vertexPropagator( [ [particle1,direction1], [particle2, direction2] ],  [newParticle, -1], propagator  ) #New particle is outgoing of this vertex
								if newPropagator != "Err" and newPropagator != None: #and !+ 0 (cos who cares?)
									# Create new list of lines to be passed in	
									newLineList = []
									newLineList.append(newLine)
									for i, line in enumerate(lineList):
										if i != index1 and i != index2:	
											newLineList.append(line)
									# Go deeper into recursion
									self.particleCalcTime += time.time() - particleStartTime
									if saveVertices:
										newVertices = vertices + [[line1, line2, newLine.getReversed()]]
									else:
										newVertices = None
									subDiagrams = self.buildDiagrams(newLineList, False, newPropagator, newPropList, newRawPropList, newVertices, fourVertices, saveVertices, optimize=optimize)
									if subDiagrams == False:
										return False
								elif newPropagator == "Err":
									print "Err: Aborting, infinite propagator created"
									self.matrixElement = 0
									self.diagramList = []
									self.diagramSetLen = 0
									return False
						# 4-vertices
						if len(lineList) > 4 and self.model.fourVertices == True:
							for index3, line3 in enumerate(lineList):
								if index3 > index2:
									newCacheCode = line1.cacheCode + line2.cacheCode + line3.cacheCode
									newCacheCodeNormed = self.getNormalizedCacheCode(newCacheCode)
									newPropList = propagators+[newCacheCodeNormed]
									newPropList.sort()
									if optimize == False:
										newRawPropList = None
									else:
										newRawPropList = rawPropagators + [newCacheCode]
									if optimize == False:
										nVert = len(newPropList)
										setLen = len(self.diagramSets)
										self.diagramSets.add(tuple(newPropList+[fourVertices+1]))
										fourVert = True
										if self.faster and len(self.diagramSets) == setLen:
											#print "False"
											fourVert = False
									else:
										#print newRawPropList
										if tuple(newRawPropList) in optimize:
											fourVert = True
										else:
											#print "False", newRawPropList
											fourVert = False
									particle3 = line3.particle
									if (particle1.particleType in particle3.couplesWith) and (particle2.particleType in particle3.couplesWith) and fourVert == True:
										if (particle3.particleType in particle1.couplesWith) and (particle3.particleType in particle2.couplesWith):
											direction3 = line3.direction
											# Try and get the line from the cache
											particleStartTime = time.time()
											placeholderParticle = particlePhysics.MasslessScalar(onShell=False)
											placeholderLine = Line( placeholderParticle, 1 ) 
											placeholderLine.cacheCode = newCacheCode
											cacheLine = self.cache[newCacheCode]
											# If it's not there, calculate it and add it
											if cacheLine == None:		
												placeholderParticle.threeMomentum.x = (particle1.threeMomentum.x*direction1 + particle2.threeMomentum.x*direction2 + particle3.threeMomentum.x*direction3) 
												placeholderParticle.threeMomentum.y = (particle1.threeMomentum.y*direction1 + particle2.threeMomentum.y*direction2 + particle3.threeMomentum.y*direction3) 
												placeholderParticle.threeMomentum.z = (particle1.threeMomentum.z*direction1 + particle2.threeMomentum.z*direction2 + particle3.threeMomentum.z*direction3) 
												placeholderParticle.energy = (particle1.energy * direction1 + particle2.energy * direction2 + particle3.energy * direction3)
												self.addToCache(placeholderLine)
											# If it is, copy it.
											else:
												placeholderParticle.threeMomentum = cacheLine.particle.getThreeMomentum()
												placeholderParticle.energy = cacheLine.particle.energy
											# Get the 4th particle type
											listNewLines = []
											vertexCharge = particle1.charge*line1.direction+particle2.charge*line2.direction+particle3.charge*line3.direction
											for particleType in particle1.couplesWith.intersection(particle2.couplesWith).intersection(particle3.couplesWith):	
												if particleType != particlePhysics.Phi and ( (particle1.particleType == particleType) and (particle2.particleType == particleType) and (particle3.particleType == particleType) ):
													self.particleCalcTime += time.time() - particleStartTime
													continue
												# Creating new particle, and a new incoming line for the rest of the diagram
												newParticle = particleType(onShell=False) # This is where I need to be checking the momentum cache instead of calculating
												newParticle.threeMomentum = placeholderParticle.threeMomentum
												newParticle.energy = placeholderParticle.energy
												# Check conservation laws are obeyed
												if newParticle.charge != vertexCharge:
													self.particleCalcTime += time.time() - particleStartTime
													continue	
												if (PRINT_ALL): print newLineLists
												newLine = Line( newParticle, 1)
												newLine.cacheCode = newCacheCode
												listNewLines.append(newLine)
											for newLine in listNewLines:
												newParticle = newLine.particle	
												newPropagator = self.model.vertexPropagator( [ [particle1,direction1], [particle2, direction2], [particle3, direction3]],  [newParticle, -1], propagator ) #New particle is outgoing of this vertex
												if newPropagator != False and newPropagator != None:
													newLineList = []
													for i, line in enumerate(lineList):
														if i != index1 and i != index2 and i != index3:	
															newLineList.append(line)
													newLineList.append(newLine)
													# Go deeper into recursion
													self.particleCalcTime += time.time() - particleStartTime
													if saveVertices:
														newVertices = vertices + [[line1, line2, line3, newLine.getReversed()]]
													else:
														newVertices = None
													subDiagrams = self.buildDiagrams(newLineList, False, newPropagator, newPropList, newRawPropList, newVertices, fourVertices+1, saveVertices, optimize=optimize)
													if subDiagrams == False:
														return False
												elif newPropagator == "Err":
													print "Err: Aborting, infinite propagator created"
													self.matrixElement = 0
													self.diagramList = []
													self.diagramSetLen = 0
													return False
			if PRINT and first: 
				print "Total time: %s seconds" %(time.time() - startTime)
				print "Internal particle time: %s seconds" %(self.particleCalcTime)
				print "Propagating time: %s seconds" %(self.diagCalcTime)
				print "Sorting time: %s seconds" %(self.sortingTime)
				print "TIMER %s " %particlePhysics.TIMER
			return
		return
		


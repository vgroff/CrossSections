from __future__ import division
from particlePhysics import Model, Particle, FourMomentum, ThreeVector
import numpy as np
import math
import random
import itertools
import traceback

class SED(Model):
	def __init__(self, samplePols=False, scramble=False, fake=False):
		Model.__init__(self, fourVertices=True)
		self.samplePols = samplePols
		self.polList    = None
		self.photons    = None
		self.iteration  = 0
		self.nPhotOut   = None
		self.events = 0
		self.MSq = 0
		self.MSqSq = 0
		self.nEvents = 0
		self.summedOver = False
		self.fake = fake
		self.scramble = scramble
		
	def vertexPropagator(self, particles, newParticle, currentPropagator, trace=False):
		electronCharge = 1
		positiveLeptons = []
		negativeLeptons = []
		leptonTypes = [ScalarElectron, ScalarPositron, ScalarMuon, ScalarAntiMuon]
		leptons = []
		#nLeptons = 0
		photons  = []
		incoming = []
		outgoing = []
		internal = []
		#print self.fake
		for particle in particles + [newParticle]:
			if particle == None:
				continue
			if particle[0].onShell:
				if particle[0] == 1:
					outgoing.append(particle)
				else:
					incoming.append(particle)
			else:
				internal.append(particle)
			particleType = particle[0].particleType
			if particleType in leptonTypes:
				if particle[0].charge == -1:
					negativeLeptons.append(particle)
					leptons.append(particle)
				elif particle[0].charge == 1:
					positiveLeptons.append(particle)
					leptons.append(particle)
			elif particleType == Photon:
				photons.append(particle)
			else:
				return "Err"
		if newParticle != None and self.fake == False:
			internalProp = newParticle[0].getPropagator(False)
			if internalProp == None:
				return "Err"
			else:
				internalProp *= 1j
		else:
			internalProp = 1
		if len(leptons) == 2 and len(photons) == 1:
			threeMomentum = ThreeVector(0,0,0)
			if len(negativeLeptons) == 2:
				energy = negativeLeptons[0][0].energy + negativeLeptons[1][0].energy
				threeMomentum.x = ( negativeLeptons[0][0].threeMomentum.x + negativeLeptons[1][0].threeMomentum.x )
				threeMomentum.y = ( negativeLeptons[0][0].threeMomentum.y + negativeLeptons[1][0].threeMomentum.y )
				threeMomentum.z = ( negativeLeptons[0][0].threeMomentum.z + negativeLeptons[1][0].threeMomentum.z )
			elif len(positiveLeptons) == 2:
				energy = -1 * (positiveLeptons[0][0].energy + positiveLeptons[1][0].energy)
				threeMomentum.x = -1*( positiveLeptons[0][0].threeMomentum.x + positiveLeptons[1][0].threeMomentum.x )
				threeMomentum.y = -1*( positiveLeptons[0][0].threeMomentum.y + positiveLeptons[1][0].threeMomentum.y )
				threeMomentum.z = -1*( positiveLeptons[0][0].threeMomentum.z + positiveLeptons[1][0].threeMomentum.z )			
			else:
				energy = negativeLeptons[0][0].energy - positiveLeptons[0][0].energy
				threeMomentum.x = ( negativeLeptons[0][0].threeMomentum.x - positiveLeptons[0][0].threeMomentum.x )
				threeMomentum.y = ( negativeLeptons[0][0].threeMomentum.y - positiveLeptons[0][0].threeMomentum.y )
				threeMomentum.z = ( negativeLeptons[0][0].threeMomentum.z - positiveLeptons[0][0].threeMomentum.z )		
			if newParticle != None and newParticle[0].particleType == Photon:
				energy *= -1j*electronCharge
				threeMomentum.scale(-1j*electronCharge)
				newParticle[0].previousPropagatorMom = [threeMomentum]
				newParticle[0].previousPropagatorE =  [energy]
				if (trace): print "Curr+Int:", currentPropagator, internalProp, "prop", threeMomentum, energy
				return currentPropagator*internalProp # Doing this vertex later!
			photon = photons[0][0]
			if photon.onShell:
				polarisation = photon.getPolarisation()[0]
				return -1j*electronCharge*(-1*threeMomentum.dot(polarisation) ) * currentPropagator * internalProp # Top row of polarisation always 0 anyway so only need -p^2
			else:
				if photon.previousPropagatorE != None:
					propSum = 0
					for i in range(len(photon.previousPropagatorE)):
						propSum += -1j*electronCharge*(photon.previousPropagatorE[i]*energy - photon.previousPropagatorMom[i].dot(threeMomentum) )
					if (trace): print "Curr+Int:", currentPropagator, internalProp, "3mom + 3mom", threeMomentum, photon.previousPropagatorMom, propSum
					return propSum * currentPropagator * internalProp
				else:
					return None
		elif len(leptons) == 2 and len(photons) == 2:
			# 4-vertex
			fourVertexCoeff = 2j*electronCharge**2
			photon1 = photons[0][0]	
			photon2 = photons[1][0]
			if (photon1.onShell) and (photon2.onShell):
				photon1Pols = photon1.getPolarisation()
				photon2Pols = photon2.getPolarisation()
				prop = 0
				for pol1 in photon1Pols:
					for pol2 in photon2Pols:
						prop += ( -1 * pol1.dot(pol2) )
				return fourVertexCoeff * prop * internalProp * currentPropagator
			else:
				if (photon1.onShell) or (photon2.onShell):
					if photon1.onShell:
						firstPhoton = photon1
						secondPhoton = photon2
					else:
						firstPhoton = photon2
						secondPhoton = photon1					
					photon1Pols = firstPhoton.getPolarisation()
					if secondPhoton.previousPropagatorE != None:
						propSum = 0
						for i in range(len(secondPhoton.previousPropagatorE)):
							for pol in photon1Pols:
								propSum += (- secondPhoton.previousPropagatorMom[i].dot(pol) )
						return  fourVertexCoeff * propSum * internalProp * currentPropagator 
					else:
						secondPhoton.previousPropagatorE = [0]
						secondPhoton.previousPropagatorMom = [pol.getScaled(fourVertexCoeff) for pol in photon1Pols] 
						if newParticle == None: return None
						else: return internalProp * currentPropagator # Doing this vertex later
				elif photon1.previousPropagatorE != None and photon2.previousPropagatorE != None:
					propSum = 0
					for i in range(len(photon1.previousPropagatorE)):
						for j in range(len(photon2.previousPropagatorE)):
							propSum += (photon1.previousPropagatorE[i]*photon2.previousPropagatorE[j] - photon1.previousPropagatorMom[i].dot(photon2.previousPropagatorMom[j]))
							if trace: print "4-v", internalProp, photon1.previousPropagatorE[i], photon2.previousPropagatorE[j], photon1.previousPropagatorMom[i], photon2.previousPropagatorMom[j]
					return fourVertexCoeff * internalProp * propSum * currentPropagator
				elif photon1.previousPropagatorE != None or photon2.previousPropagatorE != None:
					if photon1.previousPropagatorE != None:
						photonCarried = photon1
						photonNew = photon2
					else:
						photonCarried = photon2
						photonNew = photon1
					photonNew.previousPropagatorMom = [mom.getScaled(fourVertexCoeff) for mom in photonCarried.previousPropagatorMom]
					photonNew.previousPropagatorE =  [e * fourVertexCoeff for e in photonCarried.previousPropagatorE]
					return internalProp * currentPropagator					
				else:
					return None
		return None

	def getExtraCoefficients(self, incoming, outgoing):
		# Average over polarisations
		polarisationCoeff = 1
		for particle in incoming:
			if particle.particleType == Photon: 
				polarisationCoeff *= 0.5
		return polarisationCoeff
	
	def sumOver(self, incoming, outgoing):
		if self.samplePols:
			if self.polarisationN2 > 0:
				ans = self.samplePolarisations(incoming, outgoing) # Set the polarisations
				self.polarisationN2 += 1
				return ans
			else:
				self.polarisationN2 += 1 # Polarisations already set for i=0  by onNewInteraction calling setupSamplePols with self.polarisationN2 = None
				return True
		else:
			return self.sumPols(incoming, outgoing)
			
	def setupSamplePols(self, incoming, outgoing):
		if self.photons == None:
			self.polarisationN = None
			self.photons = []
			for particle in incoming+outgoing:
				if particle.particleType == Photon:
					self.photons.append(particle)
			self.polList = []
			self.polChoices([False for i in range(len(self.photons))])
			polList = self.polList
			self.nPols = len(polList)
			polN = [[] for i in range(len(self.photons)+1)]
			for pol in polList:
				polN[pol.count(True)].append(pol)
			self.polList = [ { "pols": polN[j], "prob":1/(len(self.photons)+1), "iteration":0, "sum":0, "sumSq":0, "sumR":0, "sumSqR":0 } for j in range(len(self.photons)+1)]
		if self.polarisationN2 == None:
			self.polarisationN2 = 0
			self.iteration += 1
			startPoint = int(30000 / 8)
			if self.scramble and (self.iteration == (self.nPols * startPoint) or ( self.iteration > (self.nPols * startPoint) and self.iteration % (500) == 0)):
				total = 0
				for pol in self.polList:
					var = pol["sumSq"]/pol["iteration"] - (pol["sum"]/pol["iteration"])**2 
					newProb = (var/(pol["iteration"]-1))**0.5# * 1/(len(self.photons)+1)	
					total += newProb
					pol["prob"] = newProb
					pol["error"] = newProb
				for pol in self.polList:
					pol["prob"] = pol["prob"] / total
					print pol
			choice  = random.random()
			probSum = 0
			for index, pol in enumerate(self.polList):
				prob = pol["prob"]
				if choice <= prob + probSum:
					self.polarisationN  = index
					for i in range(len(self.photons)):
						self.photons[i].setPolarisation( pol["pols"][0][i] ) 
					break
				else: 
					probSum += prob
			self.polList[self.polarisationN]["iteration"] += len(self.polList[self.polarisationN]["pols"])
			self.polList[self.polarisationN]["tempMSq"] = 0
			return True
	
	def samplePolarisations(self, incoming, outgoing):
		if self.polarisationN2 < len(self.polList[self.polarisationN]["pols"]):
			for i in range(len(self.photons)):
				self.photons[i].setPolarisation( self.polList[self.polarisationN]["pols"][self.polarisationN2][i] ) 
			return True
		else:
			self.polList[self.polarisationN]["sum"] += self.polList[self.polarisationN]["tempMSq"]
			self.polList[self.polarisationN]["sumSq"] += self.polList[self.polarisationN]["tempMSq"]**2
			self.polList[self.polarisationN]["tempMSq"] = 0
			return False
		
	def diagramCalculated(self, incoming, outgoing, diagramMSq, weight):
		if self.samplePols:
			self.polList[self.polarisationN]["tempMSq"] += diagramMSq * weight
		else:
			self.tempMSq += diagramMSq * weight
			self.nEvents += 1
			index = self.polList[self.polarisationN-1].count(True)
			self.polStats[index]["sum"] += diagramMSq * weight
			self.polStats[index]["mSq"] = self.polStats[index]["sum"] / self.polStats[index]["iteration"]
			
	def getMSq(self):
		if self.samplePols:
			MSq = 0
			n = len(self.polList) - 1
			errorSq = 0
			error2Sq = 0
			for index, pol in enumerate(self.polList):
				if pol["iteration"] > 0:
					iterations = pol["iteration"]
					mSq = pol["sum"]/iterations
					MSq += len(pol["pols"])*mSq
					errorRaw = ( (pol["sumSq"]/iterations - (pol["sum"]/iterations)**2) / (iterations - 1) )**0.5
					errorSq += (errorRaw * len(pol["pols"]))**2
			return MSq, errorSq**0.5
		else:
			#print self.MSqSq / self.events - (self.MSq / self.events)**2
			sum1 = 0
			sum2 = 0
			sumEv = 0
			n = len(self.photons)
			return self.MSq / self.events, ( (self.MSqSq / self.events - (self.MSq / self.events)**2) / (self.events - 1))**0.5


	def getNEvents(self):
		if self.samplePols:
			n = 0
			for index,pol in enumerate(self.polList):
				n += pol["iteration"]
			return n
		else:
			return self.nEvents

	def onNewInteraction(self, incoming, outgoing):
		self.setupCacheSystem(incoming, outgoing)
		if self.samplePols: 
			self.polarisationN2 = None
			self.summedOver = False
			self.setupSamplePols(incoming, outgoing)
		else: 
			self.events += 1
			self.polarisationN = 0
			self.tempMSq = 0
			self.setupSumPols(incoming, outgoing)
			for pol in self.polStats: pol["iteration"] += len(pol["pols"])
		return
		
	def setupSumPols(self, incoming, outgoing):
		if self.photons == None:
			self.photons = []
			self.polList = []
			for particle in incoming+outgoing:
				if particle.particleType == Photon:
					self.photons.append(particle)
			self.polChoices([False for i in range(len(self.photons))])
			polN = [[] for i in range(len(self.photons)+1)]
			for pol in self.polList:
				polN[pol.count(True)].append(pol)
			self.polStats = [{"sum": 0, "sumSq":0,"iteration":0, "mSq":0, "pols":polN[i]} for i in range(len(self.photons)+1)]
			
	def sumPols(self, incoming, outgoing):
		if len(self.polList) > self.polarisationN:
			for index, pol in enumerate(self.polList[self.polarisationN]):
				self.photons[index].setPolarisation(pol)
			self.polarisationN += 1
			return True
		else:
			self.MSq += self.tempMSq # / len(self.polList)
			self.MSqSq += self.tempMSq**2
			return False
		
	def polChoices(self, polList, n=0): 
		if len(polList) == n:
			self.polList.append(polList)
		else:
			polList[n] = True
			self.polChoices( list(polList), n+1)
			polList[n] = False
			self.polChoices( list(polList), n+1)
		
			

class ScalarElectron(Particle):
	i = 0
	def __init__(self, onShell = True):
		couplesWith = set([ScalarPositron, ScalarElectron, Photon])
		self.particleType = ScalarElectron
		self.antiParticleType = ScalarPositron
		self.particleType.i += 1
		self.ID = self.particleType.i
		self.model = SED
		Particle.__init__(self, 0, -1, couplesWith, onShell=onShell, name="Electron"+str(self.particleType.i))	
	

		
class ScalarPositron(Particle):
	i = 0
	def __init__(self, onShell = True):
		couplesWith = set([ScalarPositron, ScalarElectron, Photon])
		self.particleType = ScalarPositron
		self.antiParticleType = ScalarElectron
		self.particleType.i += 1
		self.ID = self.particleType.i
		self.model = SED
		Particle.__init__(self, 0, 1, couplesWith, onShell=onShell, name="Positron"+str(self.particleType.i))	
		
class ScalarMuon(Particle):
	i = 0
	def __init__(self, onShell = True):
		couplesWith = set([ScalarAntiMuon, ScalarMuon, Photon])
		self.particleType = ScalarMuon
		self.antiParticleType = ScalarAntiMuon
		self.particleType.i += 1
		self.ID = self.particleType.i
		self.model = SED
		Particle.__init__(self, 0, -1, couplesWith, onShell=onShell, name="Muon"+str(self.particleType.i))	
		
		
class ScalarAntiMuon(Particle):
	i = 0
	def __init__(self, onShell = True):
		couplesWith = set([ScalarAntiMuon, ScalarMuon, Photon])
		self.particleType = ScalarAntiMuon
		self.antiParticleType = ScalarMuon
		self.particleType.i += 1
		self.ID = self.particleType.i
		self.model = SED
		Particle.__init__(self, 0, 1, couplesWith, onShell=onShell, name="Anti-Muon"+str(self.particleType.i))	

class Photon(Particle):
	i = 0
	def __init__(self, onShell = True):
		couplesWith = set([ScalarPositron, ScalarElectron, Photon, ScalarAntiMuon, ScalarMuon])
		self.particleType = Photon
		self.antiParticleType = Photon
		self.particleType.i += 1
		self.ID = self.particleType.i
		self.model = SED
		Particle.__init__(self, 0, 0, couplesWith, onShell=onShell, name="Photon"+str(self.particleType.i))	
		self.polarisations = None
		self.previousPropagatorE = None
		self.previousPropagatorMom = None
		
	# Return polarisation set by setPolarisation/getPolarisation2
	def getPolarisation(self):
		#if self.polarisations != None: return self.polarisations[0], self.polarisations[1]
		#print "(", polarisation.x, polarisation.y, polarisation.z,")","  (", polarisation2.x, polarisation2.y, polarisation2.z
		#print polarisation.dot(polarisation2), polarisation.dot(self.threeMomentum), polarisation2.dot(self.threeMomentum)
		return [self.polarisation]
		
	# Produce 2 polarisation vectors perpendicular to direction of propagation
	def getPolarisation2(self):
		sqrt2 = 2**0.5
		polarisation  = ThreeVector(1, 0, 0)
		polarisation2 = ThreeVector(0, 1, 0)
		vector = self.threeMomentum
		#r = vector.mod()
		theta, phi = vector.getAngles()
		if vector.x != 0: direction = -vector.x/abs(vector.x)
		else: direction = -1
		polarisation.rotate( -theta, ThreeVector(0,direction*1,0) )
		polarisation.rotate( phi, ThreeVector(0,0,1) )
		polarisation2.rotate( -theta, ThreeVector(0,direction*1,0) )
		polarisation2.rotate( phi, ThreeVector(0,0,1) )
		self.polarisation = [polarisation, polarisation2]
		return polarisation, polarisation2
		
	def setPolarisation(self, up):
		sqrt2 = 2**0.5
		polarisation  = ThreeVector(1/sqrt2,  1j/sqrt2, 0)
		polarisation2 = ThreeVector(1/sqrt2, -1j/sqrt2, 0)
		#polarisation  = ThreeVector(1,0, 0)
		#polarisation2 = ThreeVector(0,1, 0)
		vector = self.threeMomentum
		#r = vector.mod()
		theta, phi = vector.getAngles()
		if vector.x != 0: direction = -vector.x/abs(vector.x)
		else: direction = -1
		#polarisation.rotate( -theta, ThreeVector(0,direction*1,0) )
		orthogVector = self.threeMomentum.cross( ThreeVector(0,0,1) )
		polarisation.rotate(-theta, orthogVector)
		#polarisation.rotate( phi, ThreeVector(0,0,1) )
		#polarisation2.rotate( -theta, ThreeVector(0,direction*1,0) 
		#~ polarisation2.rotate(theta, ThreeVector(0,1,0))
		#~ polarisation2.rotate( phi, ThreeVector(0,0,1) )
		polarisation2.rotate(-theta, orthogVector)
		self.polarisation = [polarisation, polarisation2]
		if up: self.polarisation = polarisation
		else: self.polarisation = polarisation2
		
		

from __future__ import division
import matplotlib.pyplot as plt
import copy
import numpy as np
import random
import time
import traceback

SPEED_OF_LIGHT = 3*10**8
SPEED_OF_LIGHT = 1
AU = 1.660539 * 10**-27

TIMER = 0

def getGamma(speed):
	return (1-(speed/SPEED_OF_LIGHT)**2)**(-0.5)

class ThreeVector():
	
	# Initialize the vector
	def __init__(self, x=None, y=None, z=None, r=None, theta=None, phi=None):
		if x != None and y != None and z != None:
			self.x = x
			self.y = y
			self.z = z
		elif r != None and theta != None and phi != None:
			self.x = r * np.sin(theta) * np.cos(phi)
			self.y = r * np.sin(theta) * np.sin(phi)
			self.z = r * np.cos(theta)	
		else:
			self.x = 0
			self.y = 0
			self.z = 0
	def __str__(self):
		if isinstance(self.x, complex):
			return "(%.4g %.4gj, %.4g %.4gj, %.4g %.4gj)" %(self.x.real, self.x.imag, self.y.real, self.y.imag, self.z.real, self.z.imag)
		else:
			return "(%.4g, %.4g, %.4g)" %(self.x, self.y, self.z)
	def __repr__(self):
		return self.__str__()
	# Copies the vector
	def copy(self):
		return ThreeVector(self.x, self.y, self.z)

	# Return the modulus of the vector
	def mod(self):
		return (self.x*self.x.conjugate() + self.y*self.y.conjugate() + self.z*self.z.conjugate())**0.5		
		
	# Scales the vector	
	def scale(self, scale):
		self.x = self.x * scale
		self.y = self.y * scale
		self.z = self.z * scale	
		
	# Scales the vector	
	def getScaled(self, scale):
		x = self.x * scale
		y = self.y * scale
		z = self.z * scale
		return ThreeVector(x=x, y=y, z=z)
			
	# Returns the dot product of this vector with a second vector
	def dot(self, vector2):
		return self.x * vector2.x + self.y * vector2.y + self.z * vector2.z		
			
	# Returns the cross product of this vector with a second vector
	def cross(self, vector2):
		x = self.y * vector2.z - self.z * vector2.y
		y = self.z * vector2.x - self.x * vector2.z
		z = self.x * vector2.y - self.y * vector2.x 
		return ThreeVector(x=x, y=y, z=z)
	
	# Returns the sum of two vectors	
	def add(self, vector2):
		return ThreeVector(self.x + vector2.x, self.y + vector2.y, self.z + vector2.z)
	
	# Returns the difference between this vector and another
	def subtract(self, vector2):
		return ThreeVector(self.x - vector2.x, self.y - vector2.y, self.z - vector2.z)
		
	# Returns vector to another vector (reversal of subtract)	
	def directionTo(self, vector2):
		return vector2.subtract(self)
		
	# Randomizes the vector direction, given a magnitude	
	def randomize(self, r):
		cosTheta = ( random.random() * 2 ) - 1
		theta = np.arccos(cosTheta)
		phi = random.random() * 2 * np.pi
		self.x = r * np.sin(theta) * np.cos(phi)
		self.y = r * np.sin(theta) * np.sin(phi)
		self.z = r * np.cos(theta)
		
	# Returns the spherical coordinate angles
	def getAngles(self):
		r = self.mod()
		if (r):
			theta = np.arccos( self.z/r )
			if (self.x):
				phi = np.arctan( self.y/self.x )
			else:
				if (self.y > 0):
					phi = np.pi/2
				elif (self.y < 0):
					phi = -np.pi/2
				else: 
					phi = 0
		else:
			theta = 0
			phi = 0
		return theta, phi

	
	# Rotates the vector in place
	def rotate(self, theta, axisVector):
		axisVector = axisVector.getUnit()
		x = axisVector.x
		y = axisVector.y
		z = axisVector.z
		cos = np.cos(theta)
		sin = np.sin(theta)
		rotationMatrix = [ 	ThreeVector( cos + x**2*(1-cos), x*y*(1-cos)-z*sin, x*z*(1-cos)+y*sin ),
							ThreeVector( y*x*(1-cos)+z*sin, cos+y**2*(1-cos), y*z*(1-cos)-x*sin ),
							ThreeVector( z*x*(1-cos)-y*sin, y*z*(1-cos)+x*sin, cos+z**2*(1-cos) ) ]	
		x = rotationMatrix[0].dot(self)
		y = rotationMatrix[1].dot(self)
		z = rotationMatrix[2].dot(self)	
		self.x = x
		self.y = y
		self.z = z		
		
	# Returns a new vector which is the result of this vector being rotated
	def getRotated(self, theta, axisVector):
		axisVector = axisVector.getUnit()
		x = axisVector.x
		y = axisVector.y
		z = axisVector.z
		cos = np.cos(theta)
		sin = np.sin(theta)
		rotationMatrix = [ 	ThreeVector( cos + x**2*(1-cos), x*y*(1-cos)-z*sin, x*z*(1-cos)+y*sin ),
							ThreeVector( y*x*(1-cos)+z*sin, cos+y**2*(1-cos), y*z*(1-cos)-x*sin ),
							ThreeVector( z*x*(1-cos)-y*sin, y*z*(1-cos)+x*sin, cos+z**2*(1-cos) ) ]
		return ThreeVector(rotationMatrix[0].dot(self), rotationMatrix[1].dot(self), rotationMatrix[2].dot(self))
	
	# rotate vector "onto" a new vector (i.e. apply the same rotation to it as if vector was pointing along z)
	def rotateOnto(self, vector):
		theta, phi = vector.getAngles()
		orthogVector = vector.cross( ThreeVector(0,0,1) )
		self.rotate(-theta, orthogVector)
		
			
	# Returns a new vector which is the unit vector form of this vector
	def getUnit(self):
		mod = self.mod()
		if mod != 0:
			return ThreeVector(self.x/mod, self.y/mod, self.z/mod)
		else:
			return ThreeVector(0, 0, 0)
	
	# Print the vector
	def printVector(self):
		print "(%s, %s, %s)" %(self.x, self.y, self.z)

class FourMomentum():
	#ID = 0
	def __init__(self, restMass, onShell=True):
		#Particle.ID += 1
		#self.ID = Particle.ID
		self.restMass = restMass
		self.onShell = onShell
		self.threeMomentum = ThreeVector(0,0,0)
		if (self.onShell == True):
			self.calculateEnergy()		
		
	def setMomentum(self, momentumXYZ = None, momMag=None, momTheta=None, momPhi=None):
		if momentumXYZ != None:
			self.threeMomentum = momentumXYZ.copy()
			self.threeMomentum.scale(SPEED_OF_LIGHT)
		elif momMag != None and momTheta != None and momPhi != None:
			self.threeMomentum = ThreeVector(r=momMag*SPEED_OF_LIGHT, theta=momTheta, phi=momPhi)
		if (self.onShell == True): # if restMass has been set to -1, this 4 vector does not represent an on-shell particle 
			self.calculateEnergy()
			#if (self.restMass > 0):
				#gamma = self.energy / (self.restMass * SPEED_OF_LIGHT**2)
				#conversionFactor = gamma*self.restMass*SPEED_OF_LIGHT
				#self.velocity = ThreeVector(x=self.threeMomentum.x/conversionFactor, y=self.threeMomentum.y/conversionFactor, z=self.threeMomentum.z/conversionFactor)
	
	#~ def setVelocity(self, velXYZ=None, velMag=None, velTheta=None, velPhi=None):
		#~ if velXYZ != None:
			#~ self.velocity = velXYZ.copy()
		#~ elif velMag != None and velTheta != None and velPhi != None:
			#~ self.velocity = ThreeVector(r=velMag, theta=velTheta, phi=velPhi)
		#~ if self.velocity.mod() > SPEED_OF_LIGHT:
			#~ print "ERROR HAS OCCURED, V > C"
			#~ return False
		#~ elif (self.onShell == True):
			#~ gamma = getGamma(self.velocity.mod())
			#~ self.threeMomentum = self.velocity.getScaled(gamma*self.restMass*SPEED_OF_LIGHT)
			#~ self.calculateEnergy() 
			
	def getRestFrameInfo(self):
		if self.energy == 0: traceback.print_stack()#print self.energy, self.threeMomentum
		pSq   = self.dot(self)
		beta  = self.threeMomentum.mod() / self.energy
		A = abs(pSq)**0.5
		B = self.energy
		C = self.threeMomentum.mod()
		factor = A**4-A**2*B**2+A**2*C**2
		if abs(factor) < 1e-12: factor = 0
		beta = (B*C - (factor)**0.5)/(A**2 + C**2)
		gamma = ( 1 - beta**2 )**(-0.5)
		#print pSq**0.5, B*gamma - beta*gamma*C, B, C, beta, gamma
		return pSq, beta, gamma		

	def boost(self, velocityVector):
		speed = velocityVector.mod()
		beta = speed / SPEED_OF_LIGHT
		gamma = getGamma(speed)
		direction = velocityVector.getUnit()
		x = direction.x
		y = direction.y
		z = direction.z
		boostMatrix = [	[ gamma        , -gamma*beta*x   , -gamma*beta*y   , -gamma*beta*z   ],
						[ -gamma*beta*x, 1+(gamma-1)*x**2, (gamma-1)*x*y   , (gamma-1)*x*z   ],
						[ -gamma*beta*y, (gamma-1)*y*x   , 1+(gamma-1)*y**2, (gamma-1)*y*z   ],
						[ -gamma*beta*z, (gamma-1)*z*x   , (gamma-1)*z*y   , 1+(gamma-1)*z**2] ]
		newVector = [0,0,0,0]
		for row in range(len(boostMatrix)):
			newVector[row] = (boostMatrix[row][0] * self.energy + boostMatrix[row][1] * self.threeMomentum.x + boostMatrix[row][2] * self.threeMomentum.y + boostMatrix[row][3] * self.threeMomentum.z) / SPEED_OF_LIGHT
		self.setMomentum(ThreeVector( *newVector[1:] ))
		if self.onShell == False:
			self.energy = newVector[0]
	
	def getBoosted(self, velocityVector): # this is the velocity vector of the reference frame
		speed = velocityVector.mod()
		beta = speed / SPEED_OF_LIGHT
		gamma = getGamma(speed)
		direction = velocityVector.getUnit()
		x = direction.x
		y = direction.y
		z = direction.z
		boostMatrix = [	[ gamma        , -gamma*beta*x   , -gamma*beta*y   , -gamma*beta*z   ],
						[ -gamma*beta*x, 1+(gamma-1)*x**2, (gamma-1)*x*y   , (gamma-1)*x*z   ],
						[ -gamma*beta*y, (gamma-1)*y*x   , 1+(gamma-1)*y**2, (gamma-1)*y*z   ],
						[ -gamma*beta*z, (gamma-1)*z*x   , (gamma-1)*z*y   , 1+(gamma-1)*z**2] ]
		newVector = [0,0,0,0]
		for row in range(len(boostMatrix)):
			newVector[row] = (boostMatrix[row][0] * self.energy + boostMatrix[row][1] * self.threeMomentum.x + boostMatrix[row][2] * self.threeMomentum.y + boostMatrix[row][3] * self.threeMomentum.z) / SPEED_OF_LIGHT
		newFourMom = FourMomentum(self.restMass, self.onShell)
		newFourMom.setMomentum(ThreeVector( *newVector[1:]  ))
		if self.onShell == False:
			newFourMom.energy = newVector[0]
		return newFourMom

	def calculateEnergy(self):
		self.energy = ( (self.restMass*SPEED_OF_LIGHT**2)**2 + self.threeMomentum.dot(self.threeMomentum) )**0.5
		return self.energy
		
	def getThreeMomentum(self):
		return self.threeMomentum.getScaled(1/SPEED_OF_LIGHT)

	def dot(self, momentum2):
		return self.energy*momentum2.energy - self.threeMomentum.dot(momentum2.threeMomentum)
		
	def printVector(self):
		print "%s, %s, %s, %s" %(self.energy, self.getThreeMomentum().x, self.getThreeMomentum().y, self.getThreeMomentum().z)

	def randomizeThreeMomentum(self, momMag=None):
		pass
		
	
class Particle(FourMomentum):
	def __init__(self, restMass, charge, couplesWith, onShell=True, name="Particle"):
		self.couplesWith = couplesWith
		self.charge = charge
		FourMomentum.__init__(self, restMass, onShell=onShell)
		self.name = name
		
	def __str__(self):
		return str(self.particleType) + " " + str(self.ID)

	def __repr__(self):
		return str(self.particleType) + " " + str(self.ID)

	def decay(self, particleType1, particleType2):
		if (self.restMass < particle1.restMass + particle2.restMass):
			print "ERROR, decay impossible"
			return False
		else:
			phi = random.random() * 2 * np.pi
			costheta = random.random() * 2 - 1
			theta = np.arccos(costheta)
			momentum = (SPEED_OF_LIGHT / ( 2 * self.restMass )) * ( (self.restMass**2 - particle1.restMass**2 - particle2.restMass**2)**2 - 4*(particle1.restMass*particle2.restMass)**2 )**0.5
			particle1.setMomentum( ThreeVector(r=momentum, theta=theta, phi=phi) )
			particle2.setMomentum( particle1.getThreeMomentum().getScaled(-1) )
			particle1.boost(self.velocity.getScaled(-1))
			particle2.boost(self.velocity.getScaled(-1))
			return particle1, particle2
	
	def getPropagator(self, printing = True):
		coefficient = self.dot(self) - self.restMass**2
		if coefficient != 0:
			return 1/(coefficient)
		elif printing:
			print "ERR: propagator coefficient is 0, propagator -> infinity"


class Model():
	def __init__(self, fourVertices=True):
		self.fourVertices = fourVertices
		self.summedOver = False
		self.MSq = 0
	
	def vertexPropagator(self, particles, newParticle, currentPropagator):
		return 1
		
	def getExtraCoefficients(self, incoming, outgoing):
		return 1	
		
	def setupCacheSystem(self, incoming, outgoing):
		self.cacheMax = 0
		for i in range(len(incoming)+len(outgoing)):
			self.cacheMax += 2**i
		self.caches = {}
		
	def addNewCache(self, name):
		if name not in self.caches:
			self.caches[name] = [None for i in range(self.cacheMax+1)]
			
	def getFromCache(self, name, cacheCode):
		if name in self.caches:
			return self.caches[cacheCode]
		else:
			print "Err: Cache named %s does not exist, please use addNewCache function first" %(name)
		
	def addToCache(self, obj, name, cacheCode):
		if name in self.caches:
			self.caches[name][cacheCode] = obj
		else:
			print "Err: Cache named %s does not exist, please use addNewCache function first" %(name)
		
	def summedOver(self, incoming, outgoing):
		if self.summedOver: 
			return False
		else: 
			self.summedOver = False
			return True
					
	def onNewInteraction(self, incoming, outgoing):
		self.setupCacheSystem(incoming, outgoing)
		self.summedOver = False
		return
		
	def acceptEvent(self, incoming, outgoing):
		return True
		
	def diagramCalculated(self, incoming, outgoing, diagramMSq, weight):
		self.MSq += diagramMSq * weight
		self.MSqSq += (diagramMSq * weight)**2
		return 
		
	def getMSq(self):
		return self.MSq, 0
	
class PhiModel(Model):
	def __init__(self):
		Model.__init__(self, fourVertices=True)
	
class PiPlus(Particle):
	i = 0
	def __init__(self, onShell = True):
		couplesWith = set([PiMinus, PiZero, PiPlus])
		self.particleType = PiPlus
		self.i += 1
		Particle.__init__(self, 0, 1, couplesWith, onShell=onShell, name="PiPlus"+str(self.i))
		
class PiMinus(Particle):
	i = 0
	def __init__(self, onShell = True):
		couplesWith = set([PiPlus, PiZero, PiMinus])
		self.particleType = PiMinus
		self.i += 1
		Particle.__init__(self, 0, -1, couplesWith, onShell=onShell, name="PiMinus"+str(self.i))

class PiZero(Particle):
	i = 0
	def __init__(self, onShell = True):
		couplesWith = set([PiMinus, PiPlus, PiZero])
		self.particleType = PiZero
		self.i += 1
		Particle.__init__(self, 0, 0, couplesWith, onShell=onShell, name="PiZero"+str(self.i))
		
class Phi(Particle):
	i = 0
	def __init__(self, onShell = True):
		couplesWith = set([Phi])
		self.particleType = Phi
		self.antiParticleType = Phi
		self.model = Model
		Phi.i += 1
		Particle.__init__(self, 0, 0, couplesWith, onShell=onShell, name="Phi"+str(Phi.i))	
		
class Xi(Particle):
	i = 0
	def __init__(self, onShell = True):
		couplesWith = set([Phi, Xi])
		self.particleType = Xi
		self.vertexPropagator = fakeFunction2
		Phi.i += 1
		Particle.__init__(self, 1, 0, couplesWith, onShell=onShell, name="Phi"+str(Phi.i))	
		
class MasslessScalar(Particle):
	def __init__(self, onShell = True):
		self.antiParticleType = MasslessScalar
		self.particleType = MasslessScalar
		Particle.__init__(self, 0, 0, [], onShell=onShell, name="MasslessScalar")	
		
def fakeFunction(particles, newParticle, currentPropagator):
	return 1	
	
def fakeFunction2(particles, newParticle, currentPropagator):
	return 1



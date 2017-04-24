from __future__ import division
import particlePhysics
import numpy as np
import random

# takes as input the number of in and out particles, and momentum of out particles

class Rambo():
	
	def __init__(self, incomingParticles, outgoingParticles):
		self.incomingParticles = incomingParticles
		self.massive = False
		self.outgoingParticles = outgoingParticles
		self.incomingMomentumSum = particlePhysics.FourMomentum(0, onShell = False)
		self.incomingMomentumSum.energy = 0
		for i in range(len(self.incomingParticles)):
			self.incomingMomentumSum.energy += self.incomingParticles[i].energy
			self.incomingMomentumSum.setMomentum( self.incomingMomentumSum.getThreeMomentum().add(self.incomingParticles[i].getThreeMomentum()) )	
		for particle in self.outgoingParticles:
			if particle.restMass != 0:
				self.massive = True	
		pi2log = np.log(np.pi/2)
		Z = [0 for i in range(len(self.outgoingParticles)+1)]
		Z[2] = pi2log
		for i in range(3, len(self.outgoingParticles)+1):
			Z[i] = Z[i-1]+pi2log-2.*np.log(i-2);
		for i in range(3, len(self.outgoingParticles)+1):
			Z[i] = Z[i] - np.log(i-1)
		self.Z_N = Z[len(self.outgoingParticles)]
						
	def generatePoint(self):
		outgoingMomentumSum = particlePhysics.FourMomentum(0, onShell = False)
		outgoingMomentumSum.energy = 0
		for i in range(len(self.outgoingParticles)):
			c = 2 * random.random() - 1			# cos(theta)
			f = 2 * np.pi * random.random()		# phi
			q = - np.log(random.random()*random.random())
			particle = self.outgoingParticles[i]
			particle.setMomentum(momMag=q, momTheta=np.arccos(c), momPhi=f)
			outgoingMomentumSum.energy += particle.energy
			outgoingMomentumSum.setMomentum( outgoingMomentumSum.getThreeMomentum().add(particle.getThreeMomentum()) )
		w = ( self.incomingMomentumSum.dot(self.incomingMomentumSum) )**0.5
		M = ( outgoingMomentumSum.dot(outgoingMomentumSum) )**0.5
		b = outgoingMomentumSum.getThreeMomentum().getScaled(-1/M)
		gamma = outgoingMomentumSum.energy / M 
		a = 1 / (1 + gamma)
		x = w / M
		for i in range(len(self.outgoingParticles)):
			energy = self.outgoingParticles[i].energy
			bq = b.dot( self.outgoingParticles[i].getThreeMomentum() )
			self.outgoingParticles[i].setMomentum( self.outgoingParticles[i].getThreeMomentum().add( b.getScaled(energy+a*bq) ).getScaled(x) )
		if (self.massive):
			self.stabilize(self.outgoingParticles, w)
		return self.outgoingParticles
		
	def getWeight(self):
		w = ( self.incomingMomentumSum.dot(self.incomingMomentumSum) )**0.5
		weight = 1
		if self.massive:
			pass
		nOut = len(self.outgoingParticles)
		gammaNminus1 = 1
		for i in range(1, nOut-1):
			gammaNminus1 *= i
		weight *= (np.pi/2)**(nOut-1)*w**(2*nOut-4) /  ( gammaNminus1**2 * (nOut-1) * ( 2*np.pi )**(3*nOut - 4) )
		#weight *= np.exp( (2*len(outgoingParticles)-4)*np.log(w) ) / ( 2*np.pi )**(3*len(outgoingParticles) - 4)
		return weight
		
		
	def massiveWeight(w):
		maxIter = 10
		iteration = 0
		accuracy = E * 1e-14
		p2 = []
		energies = []
		massesSq = []
		totalMass = 0
		totalE = 0
		outgoingParticles = self.outgoingParticles
		for p in outgoingParticles:
			p2.append( p.getThreeMomentum().dot(p.getThreeMomentum()) )
			totalMass += p.restMass
			massesSq.append(p.restMass**2)
			energies.append(0)
			totalE += p.energy
		x = (1-(totalMass/E)**2)**0.5
		while True:
			f0  = -E
			totalE = 0
			g0  =  0 
			xSq = x**2
			for i in range(len(self.outgoingParticles)):
				energies[i] = (massesSq[i] + xSq * p2[i])**0.5
				totalE += energies[i]
				f0 += energies[i]
				g0 += p2[i] / energies[i]
			if (abs(f0) < accuracy): 
				break
			iteration += 1
			if (iteration > maxIter): 
				break
			x -= f0/(x*g0);
			
		wt2 = 1
		wt3 = 0
		
		for p in outgoingParticles:
			v  = p.threeMomentum.mod()
			wt2 *= v/p.energy
			wt3 += v**2/p.energy
		weight = np.exp( (2*len(outgoingParticles)-3)*np.log(1/x)+np.log(wt2/wt3*w) )
			
		#~ for p in outgoingParticles:
			#~ p.setMomentum( p.getThreeMomentum().getScaled(x))

	
	def stabilize(self, outgoingParticles, E):
		maxIter = 10
		iteration = 0
		accuracy = E * 1e-14
		p2 = []
		energies = []
		massesSq = []
		totalMass = 0
		totalE = 0
		for p in outgoingParticles:
			p2.append( p.getThreeMomentum().dot(p.getThreeMomentum()) )
			totalMass += p.restMass
			massesSq.append(p.restMass**2)
			energies.append(0)
			totalE += p.energy
		x = (1-(totalMass/E)**2)**0.5
		while True:
			f0  = -E
			totalE = 0
			g0  =  0 
			xSq = x**2
			for i in range(len(self.outgoingParticles)):
				energies[i] = (massesSq[i] + xSq * p2[i])**0.5
				totalE += energies[i]
				f0 += energies[i]
				g0 += p2[i] / energies[i]
			if (abs(f0) < accuracy): 
				break
			iteration += 1
			if (iteration > maxIter): 
				break
			x -= f0/(x*g0);
			
		for p in outgoingParticles:
			p.setMomentum( p.getThreeMomentum().getScaled(x))



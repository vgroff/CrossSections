from __future__ import division
from decimal import Decimal 
import math
import matplotlib.pyplot as plt
import sys

import feynmanDiagrams
import rambo
import particlePhysics
import numpy as np
import time

saveVals = []

def monteCarlo(incoming, outgoing, n, threshold = 100):
	sumMSquared = 0
	sumMSquaredSquared = 0
	generator = rambo.Rambo( incoming, outgoing )
	energy = 0
	momentum = particlePhysics.ThreeVector(0,0,0)
	types = []
	numberOf = []
	for particle in outgoing:
		if particle.particleType not in types:
			types.append(particle.particleType)
			numberOf.append(1)
		else:
			numberOf[types.index(particle.particleType)] += 1
	symmetryCoeff = 1
	# Calculate S, due to identical particles
	for number in numberOf:
		for i in range(1, number+1):
			symmetryCoeff *= 1/i
	# Calculating incoming p^2 for a reference point, divide by 100 for the threshold p^2 (arbitrary)
	for particle in incoming:
		energy += particle.energy
		momentum.x += particle.threeMomentum.x
		momentum.y += particle.threeMomentum.y
		momentum.z += particle.threeMomentum.z
	if threshold != None:
		thresholdProp = ( (energy)**2 - momentum.dot(momentum) ) / threshold
		thresholdE = energy / threshold
	model = incoming[0].model()
	outgoing = generator.generatePoint()
	model.onNewInteraction(incoming, outgoing)
	model.sumOver(incoming, outgoing);
	interaction = feynmanDiagrams.Interaction(incoming, outgoing, model)
	interaction.calcDiagrams(saveVertices=True)
	diags = set([])
	for diagram in interaction.diagrams:
		#print "\n\nDIAG"
		#diagram.traverse()
		for subDiag in diagram.produceList():
			diags.add(subDiag)
	nDiagrams = len(interaction.diagrams)
	j=0
	# Generate n points
	step = math.ceil(n / 50)
	for i in range(n):
		if i == 0: print "0 %"
		elif i % step == 0: 
			sys.stdout.flush()
			print "\r %s %%" %(i*100/n)
			if i > 3000:
				crossSec, error = model.getMSq()
				#print error
				nEvents = model.getNEvents()
				nEvents = i
				saveVals.append( [nEvents, error] )
				#print model.polStats
		outgoing = generator.generatePoint()
		allMomenta = incoming + outgoing
		model.onNewInteraction(incoming, outgoing)
		if ( threshold == None or checkMomentaValid(allMomenta, thresholdProp, thresholdE, len(incoming)) ) and model.acceptEvent(incoming, outgoing):
			interaction = feynmanDiagrams.Interaction(incoming, outgoing, model)
			mSq = 0
			summing = model.sumOver(incoming, outgoing)
			j += 1
			while summing == True:
				interaction.changed()
				diagrams = interaction.calcDiagrams(optimize=diags, diagNumber=nDiagrams)
				diagMSq = interaction.matrixElement*interaction.matrixElement.conjugate()
				model.diagramCalculated(incoming, outgoing, diagMSq, generator.getWeight()) 
				summing = model.sumOver(incoming, outgoing)
	MSquared, errorEstimate = model.getMSq()
	fluxCoeff = 4 * ( (incoming[0].dot(incoming[1]))**2 - (incoming[0].restMass*incoming[1].restMass)**2 )**0.5 #4*sqrt( (p1.p2)^2 - (m1m2)^2)
	modelCoeff = model.getExtraCoefficients(incoming, outgoing)
	crossSection = MSquared * symmetryCoeff * modelCoeff / fluxCoeff
	errorEstimate = errorEstimate * symmetryCoeff * modelCoeff / fluxCoeff
	#print symmetryCoeff, modelCoeff, fluxCoeff, MSquared, MSquared / generator.getWeight(), generator.getWeight()
	#print generator.getWeight(),  symmetryCoeff * modelCoeff / fluxCoeff
	print "M", MSquared, MSquared / generator.getWeight(), gevToPicobarns(MSquared * symmetryCoeff * modelCoeff / (fluxCoeff))
	print "sigma", crossSection
	print "error", errorEstimate
	print "In pb: %.3E +- %.1E" %( Decimal(gevToPicobarns(crossSection.real)), Decimal(gevToPicobarns(errorEstimate.real)) )
	return crossSection.real, errorEstimate.real

def checkMomentaValid(allMomenta, thresholdProp, thresholdE, nIncoming):
	# Check that all (pi+pj)^2>threshold
	for index1, mom1 in enumerate(allMomenta):
		for index2, mom2 in enumerate(allMomenta):
			if index2 > index1:
				energy = mom1.energy + mom2.energy
				direction1 = 1
				direction2 = 1
				if index1 < nIncoming: direction1 = -1
				if index2 < nIncoming: direction2 = -1
				threeMom1 = mom1.threeMomentum
				threeMom2 = mom2.threeMomentum
				newMom = particlePhysics.ThreeVector(threeMom1.x*direction1+threeMom2.x*direction2, threeMom1.y*direction1+threeMom2.y*direction2, threeMom1.z*direction1+threeMom2.z*direction2)
				if energy < thresholdE:	
					return False
				if abs( (energy)**2 - newMom.dot(newMom) ) < thresholdProp: #Do we want abs here?
					return False
	return True

def gevToPicobarns(gev):
	return gev*0.3894*10**9

def doMC(inc, out, n=1000, threshold=None):
	inc[0].setMomentum( particlePhysics.ThreeVector(1,0,0) )
	inc[1].setMomentum( particlePhysics.ThreeVector(-1,0,0) )	
	crossSection, Err = monteCarlo(inc, out, n, threshold) 
	



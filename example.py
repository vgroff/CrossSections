import feynmanDiagrams
import rambo
import particlePhysics
import SED
import monteCarlo

# Create the particles
positron = SED.ScalarPositron()
electron = SED.ScalarElectron()
photon1  = SED.Photon()
photon2  = SED.Photon()
photon3  = SED.Photon()
# Give them momentum
electron.setMomentum( particlePhysics.ThreeVector(1,0,0) )
positron.setMomentum( particlePhysics.ThreeVector(-1,0,0) )
incoming = [electron, positron]
outgoing = [photon1, photon2, photon3]
# Generate the outgoing momenta using Rambo
generator = rambo.Rambo( incoming, outgoing )
generator.generatePoint()
# Initialise the SED model
model = electron.model()

interaction = feynmanDiagrams.Interaction(incoming, outgoing, model)
interaction.calcDiagrams(saveVertices=True)


monteCarlo.doMC( [SED.ScalarElectron(), SED.ScalarPositron()], [SED.Photon() for i in range(3)], 50000, 100)


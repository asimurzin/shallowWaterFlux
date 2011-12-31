#!/usr/bin/env python

#--------------------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#---------------------------------------------------------------------------
from Foam import ref, man


#---------------------------------------------------------------------------
def readGravitationalAcceleration( runTime, mesh ):
    ref.ext_Info() << "\nReading gravitationalProperties" << ref.nl
    
    gravitationalProperties = man.IOdictionary( man.IOobject( ref.word( "gravitationalProperties" ),
                                                              ref.fileName( runTime.constant() ),
                                                              mesh,
                                                              ref.IOobject.MUST_READ_IF_MODIFIED,
                                                              ref.IOobject.NO_WRITE ) )

    g = ref.dimensionedVector( gravitationalProperties.lookup( ref.word( "g" ) ) )
    rotating = ref.Switch( gravitationalProperties.lookup( ref.word( "rotating" ) ) )
    
    if rotating:
       Omega = ref.dimensionedVector( gravitationalProperties.lookup( ref.word( "Omega" ) ) )
    else:
       Omega = ref.dimensionedVector( ref.word( "Omega" ), -ref.dimTime, ref.vector( 0,0,0 ) )
    
    magg = g.mag()
    gHat = g / magg
    
    return gravitationalProperties, g, rotating, Omega, magg, gHat


#----------------------------------------------------------------------------
def createPhi( runTime, hU, mesh ):

    ref.ext_Info() << "Reading/calculating face flux field phi\n" << ref.nl
    
    phi = man.surfaceScalarField( man.IOobject( ref.word( "phi" ),
                                                ref.fileName( runTime.timeName() ),
                                                mesh,
                                                ref.IOobject.READ_IF_PRESENT,
                                                ref.IOobject.AUTO_WRITE ),
                                  man.linearInterpolate( hU ) & man.surfaceVectorField( mesh.Sf(), man.Deps( mesh ) ) )
    return phi
    

#----------------------------------------------------------------------------
def _createFields( runTime, mesh, Omega, gHat ):
    
    ref.ext_Info() << "Reading field h\n" << ref.nl
    h = man.volScalarField( man.IOobject( ref.word( "h" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )

    ref.ext_Info() << "Reading field h0 if present\n" << ref.nl
    h0 = man.volScalarField( man.IOobject( ref.word( "h0" ),
                                           ref.fileName( runTime.findInstance( ref.fileName( ref.word( "polyMesh" ) ), ref.word( "points" ) ) ) ,
                                           mesh,
                                           ref.IOobject.READ_IF_PRESENT ),
                             mesh,
                             ref.dimensionedScalar( ref.word( "h0" ), ref.dimLength, 0.0 ) )
    
    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ),
                            mesh )
    
    ref.ext_Info() << "Creating field hU\n" << ref.nl
    hU = man.volVectorField( man.IOobject( ref.word( "hU" ),
                                           ref.fileName( runTime.timeName() ),
                                           mesh ),
                             h * U,
                             U.ext_boundaryField().types() )
    
    ref.ext_Info() << "Creating field hTotal for post processing\n" << ref.nl
    hTotal = man.volScalarField( man.IOobject( ref.word( "hTotal" ),
                                               ref.fileName( runTime.timeName() ),
                                               mesh, 
                                               ref.IOobject.READ_IF_PRESENT,
                                               ref.IOobject.AUTO_WRITE ),
                                 h + h0 )
                             
    hTotal.write()
    
    phi = createPhi( runTime, hU, mesh )
    
    ref.ext_Info() << "Creating Coriolis Force" << ref.nl

    F = ref.dimensionedVector( ref.word( "F" ), ( ( 2.0 * Omega ) & gHat ) * gHat )
    
    return h, h0, U, hU, hTotal, phi, F


#--------------------------------------------------------------------------------------
def CourantNo( runTime, mesh, h, phi, magg ):
    CoNum = 0.0
    meanCoNum = 0.0
    waveCoNum = 0.0
    
    
    if mesh.nInternalFaces():
        phi_mag = ref.fvc.surfaceSum( phi.mag() )
        sumPhi = phi_mag.internalField() / h.internalField()
        V = mesh.V()
        CoNum = 0.5 * ( sumPhi / V.field() ).gMax() * runTime.deltaTValue()
        V_field = V.field()
        meanCoNum = 0.5 * ( sumPhi.gSum() / V_field.gSum() ) * runTime.deltaTValue()
        # Gravity wave Courant number
        waveCoNum = 0.25 * ( ref.fvc.surfaceSum( ref.fvc.interpolate( h.sqrt() ) * mesh.magSf() ).internalField() 
                             / V.field() ).gMax() * magg.sqrt().value() * runTime.deltaTValue()
        pass

    ref.ext_Info() << "Courant number mean: " << meanCoNum  << " max: " << CoNum << ref.nl
    ref.ext_Info() << "Gravity wave Courant number max: " << waveCoNum  << ref.nl
    
    return CoNum, meanCoNum, waveCoNum
    

#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    gravitationalProperties, g, rotating, Omega, magg, gHat = readGravitationalAcceleration( runTime, mesh )

    h, h0, U, hU, hTotal, phi, F = _createFields( runTime, mesh, Omega, gHat )
    
    pimple = man.pimpleControl( mesh )

    ref.ext_Info() << "\nStarting time loop\n" << ref.nl 
    
    while runTime.loop() :
        ref.ext_Info() << "\n Time = " << runTime.timeName() << ref.nl << ref.nl
        
        CourantNo( runTime, mesh, h, phi, magg )
        
        pimple.start()
        while pimple.loop():
           
           phiv = ref.surfaceScalarField( ref.word( "phiv" ), phi() / ref.fvc.interpolate( h ) ) # mixed calculations
           
           hUEqn = ref.fvm.ddt( hU ) + ref.fvm.div( phiv, hU ) 
           
           hUEqn.relax()
           
           if pimple.momentumPredictor():

              if rotating:
                  ref.solve( hUEqn + ( F ^ hU ) == -magg * h * ref.fvc.grad( h + h0 ) )
                  pass
              else:
                  ref.solve( hUEqn == -magg * h * ref.fvc.grad( h + h0 ) ) 
                  pass
              
              # Constrain the momentum to be in the geometry if 3D geometry
              if mesh.nGeometricD() == 3 :
                 hU -= ( gHat & hU ) * gHat
                 hU.correctBoundaryConditions();
                 pass
           
           for corr in range( pimple.nCorr() ): 
               hf = ref.fvc.interpolate( h )
               rUA = 1.0 / hUEqn.A()
               ghrUAf = magg * ref.fvc.interpolate( h * rUA )
               
               phih0 = ghrUAf * mesh.magSf() * ref.fvc.snGrad( h0 )
               if rotating:
                  hU << rUA * ( hUEqn .H() - ( F ^ hU ) )
                  pass
               else:
                  hU << rUA * hUEqn.H()
                  pass
               
               phi << ( ref.fvc.interpolate( hU ) & mesh.Sf() ) + ref.fvc.ddtPhiCorr( rUA, h, hU, phi )- phih0
               
               for nonOrth in range( pimple.nNonOrthCorr() + 1):
                   hEqn = ref.fvm.ddt( h ) + ref.fvc.div( phi ) - ref.fvm.laplacian( ghrUAf, h )
                   
                   hEqn.solve( mesh.solver( h.select(pimple.finalInnerIter( corr, nonOrth ) ) ) )

                   if nonOrth == pimple.nNonOrthCorr():
                      phi += hEqn.flux()
                   pass
               
               hU -= rUA * h * magg * ref.fvc.grad( h + h0 )
               
               #Constrain the momentum to be in the geometry if 3D geometry
               if mesh.nGeometricD() == 3:
                  hU -= ( gHat & hU ) * gHat
                  pass
               
               hU.correctBoundaryConditions()
               pass
           
           pimple.increment()
           pass
        
        U == hU / h
        hTotal == h + h0

        runTime.write()
        
        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << \
              "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        
        pass

    ref.ext_Info() << "End\n" << ref.nl 

    import os
    return os.EX_OK


#--------------------------------------------------------------------------------------
from Foam import FOAM_VERSION
if FOAM_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
     import sys, os
     argv = sys.argv
     os._exit( main_standalone( len( argv ), argv ) )
     pass
   pass
else:
   ref.ext_Info() << "\n\n To use this solver it is necessary to SWIG OpenFOAM-2.0.0 or higher\n"
   pass


#--------------------------------------------------------------------------------------

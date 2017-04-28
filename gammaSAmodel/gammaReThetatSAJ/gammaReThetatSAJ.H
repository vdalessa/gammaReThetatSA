/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::incompressible::RASModels::gammaReThetatSAJ

Description
    Implementation of the gamma-ReTheta model coupled with Spalart--Allmaras
    equation for predicting laminar-to-turbulent transition on arbitrary 
    geometries and grid topologies.

    For information about how to use the turbulence model, what boundary
    conditions to apply etc. please refer to the following literature:
   
    V.D'Alessandro, S. Montelpare, R.Ricci, A. Zoppi. Numerical modeling 
    of the flow past wind turbine airfoils by means of Spalart--Allmaras
    local correlation based transition model. Energy, 2017, 10.1016/j.energy.2017.04.134

    @endverbatim

SourceFiles
    gammaReThetatSAJ.C

Authors
    Valerio D'Alessandro <v.dalessandro@univpm.it; valeriodalessandro@gmail.com>
    Sergio Montelpare <s.montelpare@unich.it>
    Renato Ricci <ricci@univpm.it>
    Andrea Zoppi <a.zoppi@univpm.it; andre.zoppi@gmail.com>


\*---------------------------------------------------------------------------*/

#ifndef gammaReThetatSAJ_H
#define gammaReThetatSAJ_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class gammaReThetatSAJ Declaration
\*---------------------------------------------------------------------------*/

class gammaReThetatSAJ
:
    public RASModel
{
    // Private data
    
	    // Convergence tolerance for calculating ReThetat
	    // (not ReThetatTilda!!)
	    static const scalar tol_;
	    // Maximum iterations
	    static const int maxIter_;

        // Model coefficients
	
	// gamma-ReThetat-model coefficients
	
	    dimensionedScalar ca1_;
	    dimensionedScalar ce1_;

	    dimensionedScalar ca2_;
	    dimensionedScalar ce2_;

	    dimensionedScalar cThetat_;

	    dimensionedScalar sigmaf_;
	    dimensionedScalar sigmaThetat_;

	    dimensionedScalar s1_;


         // SA model coeffs.
            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cv2_;
            dimensionedScalar Cr1_;
            dimensionedScalar TuInf_;
            dimensionedScalar cFlength_;
            dimensionedScalar cReThetac_;
            dimensionedScalar cf1_;
            dimensionedScalar CpF_;


        wallDist y_;

        // Fields

            volScalarField gamma_;
            volScalarField ReThetatTilda_;
            volScalarField nuTilda_;
            volScalarField nut_;

    // Private member functions

        // Empirical correlations
	
	volScalarField Flength() const;
	volScalarField ReThetac() const;

	// Model functions of gamma-ReThetat-model
	
	tmp<volScalarField> Fonset() const;
	tmp<volScalarField> Fonset1() const;
	tmp<volScalarField> Fturb() const;
	tmp<volScalarField> Freattach() const;
	tmp<volScalarField> Fwake() const;
	tmp<volScalarField> FThetat() const;

	// Empirical Functions for calculating the local transition onset
	// momentum thickness Reynolds number
	void ReThetat(volScalarField& ReThetatField) const;
	scalar ReThetatEq(scalar TuInf, scalar lambda) const;

	// Intermittency modification for separation-induced transition
	tmp<volScalarField> gammaSep() const;


        tmp<volScalarField> F1(const volScalarField& CDkOmega) const;
        tmp<volScalarField> F2() const;

        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

// Funzioni protette SA
        tmp<volScalarField> chi() const;

        tmp<volScalarField> fv1(const volScalarField& chi) const;

        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

        tmp<volScalarField> fw(const volScalarField& Stilda) const;

        tmp<volScalarField> Rt () const; 

        tmp<volScalarField> G (const volScalarField& Fonset1  ) const;
// --------------------------------------------------------------------



public:

    //- Runtime type information
    TypeName("gammaReThetatSAJ");


    // Constructors

        //- Construct from components
        gammaReThetatSAJ
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
	    const word& turbulenceModelName = turbulenceModel::typeName,
	    const word& modelName = typeName
        );


    //- Destructor
    virtual ~gammaReThetatSAJ()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for gamma
        tmp<volScalarField> DgammaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DgammaEff", (nut_/sigmaf_) + nu())
            );
        }

        //- Return the effective diffusivity for ReThetatTilda
        tmp<volScalarField> DReThetatTildaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DReThetatTildaEff", sigmaThetat_*(nut_ + nu()))
            );
        }

        tmp<volScalarField> DnuTildaEff() const
        {
            return tmp<volScalarField>
            (
                 new volScalarField("DnuTildaEff", nut_ + nu())
            );
        }

             //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;


        //- Return the intermittency
        virtual tmp<volScalarField> gamma() const
        {
            return gamma_;
        }

        //- Return the local transition momentum thickness Reynolds number
        virtual tmp<volScalarField> ReThetatTilda() const
        {
            return ReThetatTilda_;
        }

	//- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

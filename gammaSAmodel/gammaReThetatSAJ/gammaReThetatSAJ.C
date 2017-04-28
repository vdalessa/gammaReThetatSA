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

\*---------------------------------------------------------------------------*/

#include "gammaReThetatSAJ.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gammaReThetatSAJ, 0);
addToRunTimeSelectionTable(RASModel, gammaReThetatSAJ, dictionary);

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

// Tolerance and maximum iteration number for calculation of ReThetat
const scalar gammaReThetatSAJ::tol_ = 1.0e-4;
const int gammaReThetatSAJ::maxIter_ = 100;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

volScalarField gammaReThetatSAJ::Flength() const
{
   	volScalarField Flength 
        (
    	    IOobject
            (
            	"Flength",
	        runTime_.timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
            ReThetatTilda_ 
	);

	//Flength = cFlength_.value() ;  
	//Flength = 40.0 ; 

        forAll(Flength, cellI)
        {
           Flength[cellI] = min( exp( 7.168-0.01173*ReThetatTilda_[cellI]) + 0.5 , 300.0  ) ;
        }
 
	return Flength ;
}

volScalarField gammaReThetatSAJ::ReThetac() const
{
   	volScalarField ReThetac 
        (
    	    IOobject
            (
            	"ReThetac",
	        runTime_.timeName(),
	        mesh_,
	        IOobject::NO_READ,
	        IOobject::NO_WRITE
	    ),
            ReThetatTilda_ 
	);

        //Medida PhD Thesis 2014 
        forAll(ReThetac, cellI)
        {
           ReThetac[cellI] = min( 0.615*ReThetatTilda_[cellI] + 61.5, ReThetatTilda_[cellI]    ) ;  
        }
	return ReThetac;
	
}

tmp<volScalarField> gammaReThetatSAJ::Fonset() const
{
	return tmp<volScalarField>
	(
	    new volScalarField
     	    (
                IOobject
       	        (
	            "Fonset",
    	            runTime_.timeName(),
	            mesh_,
	            IOobject::NO_READ,
	            IOobject::NO_WRITE
	        ),
		 max(
                      min(max( Fonset1(),pow4(Fonset1())  ) ,4.0) - max( 2.0 -pow3(Rt()/2.5)  ,0.0)   , 0.0 
                )
   	    )
	);
}

tmp<volScalarField> gammaReThetatSAJ::Fonset1() const
{
	return sqr(y_)*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))/(cf1_.value()*nu()*ReThetac());
}

tmp<volScalarField> gammaReThetatSAJ::Fturb() const
{
	return exp(-pow4(Rt()/scalar(4.0)));
}

tmp<volScalarField> gammaReThetatSAJ::Freattach() const
{
	return exp(-pow4(Rt()/scalar(20.0)));
}

tmp<volScalarField> gammaReThetatSAJ::FThetat() const
{
        volScalarField magVort = sqrt(scalar(2))*mag(skew(fvc::grad(U_)));
        magVort = max(magVort,dimensionedScalar("smallOmega",magVort.dimensions(),SMALL));
        return min
        (
            max
            (
                1.0*exp(-pow4(magSqr(U_)/(scalar(375.0)*nu()*magVort*ReThetatTilda_))),
                scalar(1.0)-sqr((ce2_*gamma_-scalar(1.0))/(ce2_-scalar(1.0)))
            ),
            scalar(1.0)
        );
}


void gammaReThetatSAJ::ReThetat(volScalarField& ReThetatField) const
{
	scalar lambda, ReThetatOld, ReThetatNew, ReThetatTol, dUds;
	volScalarField U2gradU = (sqr(U_)&&(fvc::grad(U_)));

	forAll(ReThetatField, cellI)
	{
		int iter = 0;
                scalar TuInf = TuInf_.value() ;  
  
		dUds = U2gradU[cellI]/(sqr(max(mag(U_[cellI]),SMALL)));
		// Starting value
		ReThetatNew = max(ReThetatEq(TuInf, scalar(0)),scalar(20.0));
		ReThetatTol = ReThetatNew*tol_;


		do
		{
			ReThetatOld = ReThetatNew;
			lambda = max(
		            min(
			        sqr(ReThetatOld)*nu()()[cellI]*dUds/(sqr(max(mag(U_[cellI]),SMALL))),
				scalar(0.1)
			    ),
			    scalar(-0.1)
			);
                        ReThetatNew = max(ReThetatEq(TuInf, lambda),scalar(20.0));
			if (iter++ > maxIter_)
			{
				FatalErrorIn
				(
					 "gammaReThetatSAJ::ReThetat(volScalarField& ReThetatField) const"
				)   << "Maximum number of iterations exceeded"
				    << abort(FatalError);
			}
		} while(mag(ReThetatNew-ReThetatOld) > ReThetatTol);

		ReThetatField[cellI] = ReThetatNew;
	}

}

scalar gammaReThetatSAJ::ReThetatEq(scalar TuInf, scalar lambda) const
{
        scalar FTu;
        if(TuInf > scalar(1.3))
                FTu = scalar(331.5)*pow((TuInf-scalar(0.5658)),scalar(-0.671));
        else
                FTu = scalar(1173.51)-scalar(589.428)*TuInf+scalar(0.2196)/sqr(TuInf);
        if(lambda > scalar(0))
                return FTu*(scalar(1.0)+scalar(0.275)*(scalar(1.0)-exp(scalar(-35.0)*lambda))*exp(-scalar(2.0)*TuInf));
        else
                return FTu*(scalar(1.0)+(scalar(12.986)*lambda+scalar(123.66)*sqr(lambda)+scalar(405.689)*pow3(lambda))*exp(-pow((TuInf/scalar(1.5)),scalar(1.5))));

}


tmp<volScalarField> gammaReThetatSAJ::gammaSep() const
{
	return FThetat()*min
	(
	    s1_*Freattach()*max
	    (
	        sqr(y_)*sqrt(scalar(2.0))*mag(symm(fvc::grad(U_)))/(scalar(3.235)*nu()*ReThetac())-scalar(1.0),
		scalar(0.0)
	    ),
	    scalar(2.0)
	);
}

tmp<volScalarField> gammaReThetatSAJ::G(const volScalarField& Fonset1) const
{
  scalar max_fonset = gMax ( Fonset1) ;  
  if ( max_fonset >= 1.0 ) 
  { 
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "G",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("1", dimensionSet(0, 0, 0, 0, 0), 1.0)
        )
     );
    }

  if ( max_fonset < 1.0 )   
  {   

     return pow(Fonset1, scalar(CpF_.value())) ; 
  }
}
// Spalart--Allmaras equation
//
tmp<volScalarField> gammaReThetatSAJ::chi() const
{
    return nuTilda_/nu();
}

tmp<volScalarField> gammaReThetatSAJ::fv1(const volScalarField& chi) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cv1_));
}

tmp<volScalarField> gammaReThetatSAJ::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0/pow3(scalar(1) + chi/Cv2_);
}


tmp<volScalarField> gammaReThetatSAJ::fw(const volScalarField& Stilda) const
{
    volScalarField r
    (
        min
        (
            nuTilda_
           /(
               max
               (
                   Stilda,
                   dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
               )
              *sqr(kappa_*y_)
            ),
            scalar(10.0)
        )
    );
    r.boundaryField() == 0.0;

    volScalarField r2tmp=r;

    r2tmp = -min( dimensionedScalar("zero1", r.dimensions(), 0.0)  , r) + max( min(r,10.0) , dimensionedScalar("zero1", r.dimensions(), 0.0)   );
    const volScalarField g(r2tmp + Cw2_*(pow6(r2tmp) - r2tmp));

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}

tmp<volScalarField> gammaReThetatSAJ::Rt() const
{

    volScalarField chi3 = pow(nuTilda_/nu() , 3.0) ;
    return (chi3/(chi3 + pow3(Cv1_)))*nuTilda_/nu() ; 
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gammaReThetatSAJ::gammaReThetatSAJ
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, lamTransportModel, turbulenceModelName),

    ca1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ca1",
            coeffDict_,
            2.0
        )
    ),
    ce1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce1",
            coeffDict_,
            1.0
        )
    ),
    ca2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ca2",
            coeffDict_,
            0.06
        )
    ),
    ce2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ce2",
            coeffDict_,
            50.0
        )
    ),
    cThetat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cThetat",
            coeffDict_,
            0.03
        )
    ),
    sigmaf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaf",
            coeffDict_,
            1.0
        )
    ),
    sigmaThetat_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaThetat",
            coeffDict_,
            2.0
        )
    ),
    s1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "s1",
            coeffDict_,
            2.0
        )
    ),
     sigmaNut_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNut",
            coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    Cb1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb1",
            coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb2",
            coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            coeffDict_,
            2.0
        )
    ),
      Cv1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv1",
            coeffDict_,
            7.1
        )
    ),
    Cv2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv2",
            coeffDict_,
            5.0
        )
    ),

   Cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr1",
            coeffDict_,
            0.3
        )
    ),

    TuInf_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "TuInf",
            coeffDict_,
            1e-2
        )
    ),

    cFlength_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cFlength",
            coeffDict_,
            40.0
        )
    ),

    cReThetac_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cReThetac",
            coeffDict_,
            0.62
        )
    ),

    cf1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cf1",
            coeffDict_,
            2.193
        )
    ),

    CpF_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CpF",
            coeffDict_,
            0.5
        )
    ),

    y_(mesh_),

    gamma_
    (
        IOobject
        (
            "gamma",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_
    ),
    ReThetatTilda_
    (
        IOobject
        (
            "ReThetatTilda",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_
    ),
    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    )
{

    volScalarField nuPos  = max(nuTilda_        , dimensionedScalar("0.00001", dimensionSet(0, 2, -1, 0, 0),0.00001 ) ) ;
    nut_ = nuTilda_*(1.0/ ( scalar(1.0) + pow((Cv1_*nu()/nuPos) ,3.0)   )  );
    nut_.correctBoundaryConditions();

    printCoeffs();
}


//****************************************************************************************************


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> gammaReThetatSAJ::k() const
{
/*    WarningIn("tmp<volScalarField> SpalartAllmarasQRC::k() const")
        << "Turbulence kinetic energy not defined for Spalart-Allmaras model. "
        << "Returning zero field" << endl; */

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0)
        )
    );
}

tmp<volScalarField> gammaReThetatSAJ::epsilon() const
{
    WarningIn("tmp<volScalarField> SpalartAllmarasQRC::epsilon() const")
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0)
        )
    );
}


tmp<volSymmTensorField> gammaReThetatSAJ::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
             ((2.0/3.0)*I)*k() - nut()*twoSymm(fvc::grad(U_)) 
        )
    );
}


tmp<volSymmTensorField> gammaReThetatSAJ::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> gammaReThetatSAJ::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}

tmp<fvVectorMatrix> gammaReThetatSAJ::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}

bool gammaReThetatSAJ::read()
{
    if (RASModel::read())
    {
	ca1_.readIfPresent(coeffDict());
	ce1_.readIfPresent(coeffDict());
	ca2_.readIfPresent(coeffDict());
	ce2_.readIfPresent(coeffDict());
	cThetat_.readIfPresent(coeffDict());
	sigmaf_.readIfPresent(coeffDict());
	sigmaThetat_.readIfPresent(coeffDict());
	s1_.readIfPresent(coeffDict());
	TuInf_.readIfPresent(coeffDict());
	cFlength_.readIfPresent(coeffDict());
	cReThetac_.readIfPresent(coeffDict());
        return true;
    }
    else
    {
        return false;
    }
}


void gammaReThetatSAJ::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    const volScalarField chi(this->chi());
    const volScalarField fv1(this->fv1(chi));

    volScalarField pp =  sqrt(2.0)*mag(skew(fvc::grad(U_))) + 2.0*min( dimensionedScalar("0.0", dimensionSet(0, 0, -1, 0, 0),0.0), sqrt(2.0)*mag(symm(fvc::grad(U_))) - sqrt(2.0)*mag(skew(fvc::grad(U_)))) ;

    const volScalarField Stilda
    (
         pp + fv2(chi, fv1)*nuTilda_/sqr(kappa_*y_)   
       //max (  pp + fv2(chi, fv1)*nuTilda_/sqr(kappa_*y_) , 0.3*pp  )  
    );

    volScalarField gammaEff  = max( gamma_ , gammaSep() ); 
    volScalarField S2 = magSqr(symm(fvc::grad(U_)));


    //nuTilda Equation
    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(nuTilda_)
      + fvm::div(phi_, nuTilda_)
      - (1.0/sigmaNut_)*fvm::laplacian(DnuTildaEff(), nuTilda_)
      - Cb2_/sigmaNut_*magSqr(fvc::grad(nuTilda_))
     ==
        fvm::Sp(gammaEff*Cb1_*Stilda, nuTilda_) 
       -fvm::Sp(  max( min(gamma_ , scalar(0.5)) , 1.0  ) *  Cw1_*fw(Stilda)*nuTilda_/sqr(y_), nuTilda_)
    );
    nuTildaEqn().relax();
    solve(nuTildaEqn);
    bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    nuTilda_.correctBoundaryConditions();

      volScalarField nuPos  = max(nuTilda_        , dimensionedScalar("0.00001", dimensionSet(0, 2, -1, 0, 0),0.00001 ) ) ;
      nut_ = nuTilda_*(1.0/ ( scalar(1.0) + pow((Cv1_*nu()/nuPos) ,3.0)   )  );

    nut_.correctBoundaryConditions();


    // local transition onset momentum thickness Reynolds number
    volScalarField ReThetatField
    (
        IOobject
        (
            "ReThetatField",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ReThetatTilda_
    );
    ReThetat(ReThetatField);
    
    // OUTPUT FUNCTIONS
    if(runTime_.outputTime())
    {
   
    }

    // Transition onset momentum thickness Reynolds number equation
    tmp<fvScalarMatrix> ReThetatTildaEqn
    (
        fvm::ddt(ReThetatTilda_)
      + fvm::div(phi_, ReThetatTilda_)
      - fvm::laplacian(DReThetatTildaEff(), ReThetatTilda_)
     ==
        max(cThetat_*magSqr(U_)*(scalar(1.0)-FThetat())*ReThetatField/(scalar(500.0)*nu()) , dimensionedScalar("1e-8", dimensionSet(0, 0, -1, 0, 0), 1e-8 )  )
      - fvm::Sp(cThetat_*magSqr(U_)*(scalar(1.0)-FThetat())/(scalar(500.0)*nu()), ReThetatTilda_)
    );

    ReThetatTildaEqn().relax();
    solve(ReThetatTildaEqn);

    bound(ReThetatTilda_,scalar(20));
  

    // Intermittency equation
    
    tmp<fvScalarMatrix> gammaEqn
    (
        fvm::ddt(gamma_)
      + fvm::div(phi_, gamma_)
      - fvm::laplacian(DgammaEff(), gamma_)
     ==
        max(Flength()*ca1_*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))*sqrt(Fonset()*gamma_), dimensionedScalar("1e-8", dimensionSet(0, 0, -1, 0, 0), 1e-8 )  )
      - fvm::Sp
        (
	    Flength()*ca1_*sqrt(scalar(2))*mag(symm(fvc::grad(U_)))*sqrt(Fonset()*gamma_)*ce1_,
	    gamma_
	)
      + max(ca2_*sqrt(scalar(2))*mag(skew(fvc::grad(U_)))*Fturb()*gamma_ , dimensionedScalar("1e-8", dimensionSet(0, 0, -1, 0, 0), 1e-8 )  )
      - fvm::Sp
        (
            ce2_*ca2_*sqrt(scalar(2))*mag(skew(fvc::grad(U_)))*Fturb()*gamma_,
            gamma_
        )

    ); 

    gammaEqn().relax();
    solve(gammaEqn);

    bound(gamma_,scalar(0));


}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

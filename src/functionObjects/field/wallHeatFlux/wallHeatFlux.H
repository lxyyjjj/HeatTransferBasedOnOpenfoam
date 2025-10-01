/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::functionObjects::wallHeatFlux

Group
    grpForcesFunctionObjects

Description
    Computes the wall-heat flux at selected wall patches.

    \table
      Operand        | Type              | Location
      input          | -                 | -
      output file    | dat               | postProcessing/<FO>/<time>/field
      output field   | volScalarField (only boundaryField) <!--
                     --> | postProcessing/<FO>/<time>/outField
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    wallHeatFluxFO
    {
        // Mandatory entries
        type        wallHeatFlux;
        libs        (fieldFunctionObjects);

        // Optional entries
        patches     (<patch1> ... <patchN>); // (wall1 "(wall2|wall3)");
        qr          <word>;;

        // Inherited entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                           | Type | Reqd | Deflt
      type         | Type name: wallHeatFlux               | word |  yes  | -
      libs         | Library name: fieldFunctionObjects    | word |  yes  | -
      patches    | Names of operand patches   | wordList | no | all wall patches
      qr         | Name of radiative heat flux field       | word | no | qr
    \endtable

    The inherited entries are elaborated in:
      - \link functionObject.H \endlink
      - \link writeFile.H \endlink

SourceFiles
    wallHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_wallHeatFlux_H
#define Foam_functionObjects_wallHeatFlux_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class wallHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class wallHeatFlux
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected Data

        //- Wall patches to process (optionally filtered by name)
        labelList patchIDs_;

        //- Name of radiative heat flux name
        word qrName_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the heat-flux
        void calcHeatFlux
        (
            const volScalarField& alpha,
            const volScalarField& he,
            volScalarField& wallHeatFlux
        );


public:

    //- Runtime type information
    TypeName("wallHeatFlux");


    // Constructors

        //- Construct from name, Time and dictionary
        wallHeatFlux
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        wallHeatFlux(const wallHeatFlux&) = delete;

        //- No copy assignment
        void operator=(const wallHeatFlux&) = delete;


    //- Destructor
    virtual ~wallHeatFlux() = default;


    // Member Functions

        //- Read the function-object dictionary
        virtual bool read(const dictionary& dict);

        //- Execute the function-object operations
        virtual bool execute();

        //- Write the function-object results
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

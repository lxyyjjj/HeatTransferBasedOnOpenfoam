tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

if (Y.size())
{
    mvConvection = tmp<fv::convectionScheme<scalar>>
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );
}

{
    if (reactionPtr)
    {
        reactionPtr->solve(runTime.deltaT().value());
        Qdot = reactionPtr->Qdot();
    }
    else
    {
        Qdot = dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero);
    }

    volScalarField Yt
    (
        IOobject("Yt", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("Yt", dimless, 0)
    );

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(thermo.mu(), Yi)
             ==
                fvOptions(rho, Yi)
            );

            if (reactionPtr)
            {
                YiEqn == reactionPtr->RR(i);
            }

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.clamp_min(0);
            Yt += Yi;
        }
    }

    if (Y.size())
    {
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].clamp_min(0);
    }
}

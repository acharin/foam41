/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    particleTracks

Description
    Generates a VTK file of particle tracks for cases that were computed using
    a tracked-parcel-type cloud

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "CloudTemplate.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "foamTime.H"
#include "timeSelector.H"
#include "OFstream.H"
#include "passiveParticleCloud.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"

    #include "setRootCase.H"

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fileName vtkPath(runTime.path()/"VTK");
    mkDir(vtkPath);

    Info<< "Scanning times to determine track data for cloud " << cloudName
        << nl << endl;

    labelList maxIds(Pstream::nProcs(), -1);
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        Info<< "    Reading particle positions" << endl;
        passiveParticleCloud myCloud(mesh, cloudName);
        Info<< "    Read " << returnReduce(myCloud.size(), sumOp<label>())
            << " particles" << endl;

        forAllConstIter(passiveParticleCloud, myCloud, iter)
        {
            label origId = iter().origId();
            label origProc = iter().origProc();

            if (origProc >= maxIds.size())
            {
                // Expand size
                maxIds.setSize(origProc+1, -1);
            }

            maxIds[origProc] = max(maxIds[origProc], origId);
        }
    }

    label maxNProcs = returnReduce(maxIds.size(), maxOp<label>());

    Info<< "Detected particles originating from " << maxNProcs
        << " processors." << nl << endl;

    maxIds.setSize(maxNProcs, -1);

    Pstream::listCombineGather(maxIds, maxEqOp<label>());
    Pstream::listCombineScatter(maxIds);

    labelList numIds = maxIds + 1;

    Info<< nl << "Particle statistics:" << endl;
    forAll(maxIds, procI)
    {
        Info<< "    Found " << numIds[procI] << " particles originating"
            << " from processor " << procI << endl;
    }
    Info<< nl << endl;


    // calc starting ids for particles on each processor
    List<label> startIds(numIds.size(), 0);
    for (label i = 0; i < numIds.size()-1; i++)
    {
        startIds[i+1] += startIds[i] + numIds[i];
    }
    label nParticle = startIds[startIds.size()-1] + numIds[startIds.size()-1];



    // number of tracks to generate
    label nTracks = nParticle/sampleFrequency;

    // storage for all particle tracks
    List<DynamicList<vector> > allTracks(nTracks);

    Info<< "\nGenerating " << nTracks << " particle tracks for cloud "
        << cloudName << nl << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        List<pointField> allPositions(Pstream::nProcs());
        List<labelField> allOrigIds(Pstream::nProcs());
        List<labelField> allOrigProcs(Pstream::nProcs());

        // Read particles. Will be size 0 if no particles.
        Info<< "    Reading particle positions" << endl;
        passiveParticleCloud myCloud(mesh, cloudName);

        // collect the track data on all processors that have positions
        allPositions[Pstream::myProcNo()].setSize
        (
            myCloud.size(),
            point::zero
        );
        allOrigIds[Pstream::myProcNo()].setSize(myCloud.size(), 0);
        allOrigProcs[Pstream::myProcNo()].setSize(myCloud.size(), 0);

        label i = 0;
        forAllConstIter(passiveParticleCloud, myCloud, iter)
        {
            allPositions[Pstream::myProcNo()][i] = iter().position();
            allOrigIds[Pstream::myProcNo()][i] = iter().origId();
            allOrigProcs[Pstream::myProcNo()][i] = iter().origProc();
            i++;
        }

        // collect the track data on the master processor
        Pstream::gatherList(allPositions);
        Pstream::gatherList(allOrigIds);
        Pstream::gatherList(allOrigProcs);

        Info<< "    Constructing tracks" << nl << endl;
        if (Pstream::master())
        {
            forAll(allPositions, procI)
            {
                forAll(allPositions[procI], i)
                {
                    label globalId =
                        startIds[allOrigProcs[procI][i]]
                      + allOrigIds[procI][i];

                    if (globalId % sampleFrequency == 0)
                    {
                        label trackId = globalId/sampleFrequency;
                        if (allTracks[trackId].size() < maxPositions)
                        {
                            allTracks[trackId].append
                            (
                                allPositions[procI][i]
                            );
                        }
                    }
                }
            }
        }
    }

    if (Pstream::master())
    {
        OFstream vtkTracks(vtkPath/"particleTracks.vtk");

        Info<< "\nWriting particle tracks to " << vtkTracks.name()
            << nl << endl;

        // Total number of points in tracks + 1 per track
        label nPoints = 0;
        forAll(allTracks, trackI)
        {
            nPoints += allTracks[trackI].size();
        }

        vtkTracks
            << "# vtk DataFile Version 2.0" << nl
            << "particleTracks" << nl
            << "ASCII" << nl
            << "DATASET POLYDATA" << nl
            << "POINTS " << nPoints << " float" << nl;

        // Write track points to file
        forAll(allTracks, trackI)
        {
            forAll(allTracks[trackI], i)
            {
                const vector& pt = allTracks[trackI][i];
                vtkTracks << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
            }
        }

        // write track (line) connectivity to file
        vtkTracks << "LINES " << nTracks << ' ' << nPoints+nTracks << nl;

        // Write ids of track points to file
        label globalPtI = 0;
        forAll(allTracks, trackI)
        {
            vtkTracks << allTracks[trackI].size();

            forAll(allTracks[trackI], i)
            {
                vtkTracks << ' ' << globalPtI;
                globalPtI++;
            }

            vtkTracks << nl;
        }

        Info<< "end" << endl;
    }

    return 0;
}


// ************************************************************************* //

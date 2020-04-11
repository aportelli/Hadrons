/*
 * Modules.hpp, part of Hadrons ()
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <dc-erbe1@tesseract-login1.ib0.sgi.cluster.dirac.ed.ac.uk>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Michael Marshall <michael.marshall@ed.ac.uk>
 * Author: Nils Asmussen <n.asmussen@soton.ac.uk>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Vera Guelpers <Vera.Guelpers@ed.ac.uk>
 * Author: Vera Guelpers <vmg1n14@soton.ac.uk>
 * Author: ferben <ferben@c180030.wlan.net.ed.ac.uk>
 * Author: ferben <ferben@c183212.wlan.net.ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
 * Author: ferben <ferben@localhost.localdomain>
 * Author: fionnoh <fionnoh@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#include <Hadrons/Modules/MAction/DWF.hpp>
#include <Hadrons/Modules/MAction/MobiusDWF.hpp>
#include <Hadrons/Modules/MAction/ScaledDWF.hpp>
#include <Hadrons/Modules/MAction/Wilson.hpp>
#include <Hadrons/Modules/MAction/WilsonClover.hpp>
#include <Hadrons/Modules/MAction/ZMobiusDWF.hpp>
#include <Hadrons/Modules/MContraction/A2AAslashField.hpp>
#include <Hadrons/Modules/MContraction/A2AFourQuarkContraction.hpp>
#include <Hadrons/Modules/MContraction/A2ALoop.hpp>
#include <Hadrons/Modules/MContraction/A2AMesonField.hpp>
#include <Hadrons/Modules/MContraction/Baryon.hpp>
#include <Hadrons/Modules/MContraction/DiscLoop.hpp>
#include <Hadrons/Modules/MContraction/Gamma3pt.hpp>
#include <Hadrons/Modules/MContraction/Meson.hpp>
#include <Hadrons/Modules/MContraction/SigmaToNucleonEye.hpp>
#include <Hadrons/Modules/MContraction/SigmaToNucleonNonEye.hpp>
#include <Hadrons/Modules/MContraction/WeakEye3pt.hpp>
#include <Hadrons/Modules/MContraction/WeakMesonDecayKl2.hpp>
#include <Hadrons/Modules/MContraction/WeakNonEye3pt.hpp>
#include <Hadrons/Modules/MDistil/Distil.hpp>
#include <Hadrons/Modules/MDistil/DistilPar.hpp>
#include <Hadrons/Modules/MDistil/DistilVectors.hpp>
#include <Hadrons/Modules/MDistil/LapEvec.hpp>
#include <Hadrons/Modules/MDistil/Noises.hpp>
#include <Hadrons/Modules/MDistil/PerambFromSolve.hpp>
#include <Hadrons/Modules/MDistil/Perambulator.hpp>
#include <Hadrons/Modules/MFermion/EMLepton.hpp>
#include <Hadrons/Modules/MFermion/FreeProp.hpp>
#include <Hadrons/Modules/MFermion/GaugeProp.hpp>
#include <Hadrons/Modules/MGauge/Electrify.hpp>
#include <Hadrons/Modules/MGauge/FundtoHirep.hpp>
#include <Hadrons/Modules/MGauge/GaugeFix.hpp>
#include <Hadrons/Modules/MGauge/Random.hpp>
#include <Hadrons/Modules/MGauge/StochEm.hpp>
#include <Hadrons/Modules/MGauge/StoutSmearing.hpp>
#include <Hadrons/Modules/MGauge/Unit.hpp>
#include <Hadrons/Modules/MGauge/UnitEm.hpp>
#include <Hadrons/Modules/MIO/LoadA2AMatrixDiskVector.hpp>
#include <Hadrons/Modules/MIO/LoadA2AVectors.hpp>
#include <Hadrons/Modules/MIO/LoadBinary.hpp>
#include <Hadrons/Modules/MIO/LoadCoarseEigenPack.hpp>
#include <Hadrons/Modules/MIO/LoadCosmHol.hpp>
#include <Hadrons/Modules/MIO/LoadDistilNoise.hpp>
#include <Hadrons/Modules/MIO/LoadEigenPack.hpp>
#include <Hadrons/Modules/MIO/LoadNersc.hpp>
#include <Hadrons/Modules/MIO/LoadPerambulator.hpp>
#include <Hadrons/Modules/MNPR/Amputate.hpp>
#include <Hadrons/Modules/MNPR/Bilinear.hpp>
#include <Hadrons/Modules/MNPR/FourQuark.hpp>
#include <Hadrons/Modules/MNoise/FullVolumeSpinColorDiagonal.hpp>
#include <Hadrons/Modules/MNoise/SparseSpinColorDiagonal.hpp>
#include <Hadrons/Modules/MNoise/TimeDilutedSpinColorDiagonal.hpp>
#include <Hadrons/Modules/MScalar/ChargedProp.hpp>
#include <Hadrons/Modules/MScalar/FreeProp.hpp>
#include <Hadrons/Modules/MScalar/Scalar.hpp>
#include <Hadrons/Modules/MScalarSUN/Div.hpp>
#include <Hadrons/Modules/MScalarSUN/EMT.hpp>
#include <Hadrons/Modules/MScalarSUN/Grad.hpp>
#include <Hadrons/Modules/MScalarSUN/StochFreeField.hpp>
#include <Hadrons/Modules/MScalarSUN/TrKinetic.hpp>
#include <Hadrons/Modules/MScalarSUN/TrMag.hpp>
#include <Hadrons/Modules/MScalarSUN/TrPhi.hpp>
#include <Hadrons/Modules/MScalarSUN/TransProj.hpp>
#include <Hadrons/Modules/MScalarSUN/TwoPoint.hpp>
#include <Hadrons/Modules/MScalarSUN/TwoPointNPR.hpp>
#include <Hadrons/Modules/MScalarSUN/Utils.hpp>
#include <Hadrons/Modules/MSink/Point.hpp>
#include <Hadrons/Modules/MSink/Smear.hpp>
#include <Hadrons/Modules/MSolver/A2AAslashVectors.hpp>
#include <Hadrons/Modules/MSolver/A2AVectors.hpp>
#include <Hadrons/Modules/MSolver/Guesser.hpp>
#include <Hadrons/Modules/MSolver/LocalCoherenceLanczos.hpp>
#include <Hadrons/Modules/MSolver/MixedPrecisionRBPrecCG.hpp>
#include <Hadrons/Modules/MSolver/RBPrecCG.hpp>
#include <Hadrons/Modules/MSource/Convolution.hpp>
#include <Hadrons/Modules/MSource/Gauss.hpp>
#include <Hadrons/Modules/MSource/JacobiSmear.hpp>
#include <Hadrons/Modules/MSource/Momentum.hpp>
#include <Hadrons/Modules/MSource/MomentumPhase.hpp>
#include <Hadrons/Modules/MSource/Point.hpp>
#include <Hadrons/Modules/MSource/SeqAslash.hpp>
#include <Hadrons/Modules/MSource/SeqConserved.hpp>
#include <Hadrons/Modules/MSource/SeqGamma.hpp>
#include <Hadrons/Modules/MSource/Wall.hpp>
#include <Hadrons/Modules/MSource/Z2.hpp>
#include <Hadrons/Modules/MUtilities/PrecisionCast.hpp>
#include <Hadrons/Modules/MUtilities/RandomVectors.hpp>

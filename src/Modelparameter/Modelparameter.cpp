#include "Modelparameter.hpp"
using namespace scai;
using namespace KITGPI;

/*! \brief Getter method for partitionedIn */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getPartitionedIn()
{
    return (PartitionedIn);
}

/*! \brief Getter method for partitionedOut */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getPartitionedOut()
{
    return (PartitionedOut);
}

/*! \brief Getter method for parametrisation */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getParametrisation()
{
    return (parametrisation);
}



/*! \brief Getter method for relaxation frequency */
template <typename ValueType>
ValueType KITGPI::Modelparameter::Modelparameter<ValueType>::getRelaxationFrequency() const
{
    return (relaxationFrequency);
}

/*! \brief Set relaxation frequency
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRelaxationFrequency(ValueType const setRelaxationFrequency)
{
dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
dirtyFlagAveraging = true; // If S-Wave velocity will be changed, averaging needs to be redone
relaxationFrequency=setRelaxationFrequency;
}

/*! \brief Getter method for number of relaxation mechanisms */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::getNumRelaxationMechanisms() const
{
    return (numRelaxationMechanisms);
}

/*! \brief Set number of Relaxation Mechanisms
 */ 
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setNumRelaxationMechanisms(IndexType const setNumRelaxationMechanisms)
{
dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
dirtyFlagAveraging = true; // If S-Wave velocity will be changed, averaging needs to be redone
numRelaxationMechanisms=setNumRelaxationMechanisms;
}

/*! \brief Init a single modelparameter by a constant value
 *
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param value Value which will be used to initialize the single modelparameter to a homogenoeus model
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, ValueType value)
{
    allocateModelparameter(vector, ctx, dist);

    vector = value;
}

/*! \brief Init a single modelparameter by reading a model from an external file
 *
 *  Reads a single model from an external mtx file.
 \param vector Singel modelparameter which will be initialized
 \param ctx Context
 \param dist Distribution
 \param filename Location of external file which will be read in
 \param partitionedIn Partitioned input
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::initModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist, std::string filename, IndexType partitionedIn)
{

    allocateModelparameter(vector, ctx, dist);

    readModelparameter(vector, filename, dist, partitionedIn);

    vector.redistribute(dist);
}

/*! \brief Write singe modelparameter to an external file
 *
 *  Write a single model to an external file block.
 \param vector Single modelparameter which will be written to filename
 \param filename Name of file in which modelparameter will be written
 \param partitionedOut Partitioned output
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::writeModelparameter(scai::lama::Vector<ValueType> const &vector, std::string filename, IndexType partitionedOut) const
{
    PartitionedInOut::PartitionedInOut<ValueType> partitionOut;

    switch (partitionedOut) {
    case false:
        vector.writeToFile(filename);
	HOST_PRINT(vector.getDistributionPtr()->getCommunicatorPtr(), "writing " << filename << "\n");
        break;

    case true:
        partitionOut.writeToDistributedFiles(vector, filename);
        break;

    default:
        COMMON_THROWEXCEPTION("Unexpected output option!")
        break;
    }
};

/*! \brief Read a modelparameter from file
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::readModelparameter(scai::lama::Vector<ValueType> &vector, std::string filename, scai::dmemo::DistributionPtr dist, IndexType partitionedIn)
{

    PartitionedInOut::PartitionedInOut<ValueType> partitionIn;

    switch (partitionedIn) {
    case false:
        partitionIn.readFromOneFile(vector, filename, dist);
        break;

    case true:
        partitionIn.readFromDistributedFiles(vector, filename, dist);
        break;

    default:
        COMMON_THROWEXCEPTION("Unexpected input option!")
        break;
    }
};

/*! \brief Allocate a single modelparameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::allocateModelparameter(scai::lama::Vector<ValueType> &vector, scai::hmemo::ContextPtr ctx, scai::dmemo::DistributionPtr dist)
{
    vector.setContextPtr(ctx);
    vector.allocate(dist);
};

/*! \brief Calculate a modulus from velocity
 *
 *  Calculates Modulus = pow(Velocity,2) *  Density
 *
 \param vecVelocity Velocity-Vector which will be used in the calculation
 \param vecDensity Density-Vector which will be used in the calculation
 \param vectorModulus Modulus-Vector which is calculated
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcModulusFromVelocity(scai::lama::Vector<ValueType> &vecVelocity, scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vectorModulus)
{

    vectorModulus = vecDensity;
    vectorModulus *= vecVelocity;
    vectorModulus *= vecVelocity;
};

/*! \brief Calculate velocities from a modulus
 *
 *  Calculates Velocity = sqrt( Modulu / Density )
 *
 \param vectorModulus Modulus-Vector which will be used in the calculation
 \param vecDensity Density-Vector which will be used in the calculation
 \param vecVelocity Velocity-Vector which is calculated
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcVelocityFromModulus(scai::lama::Vector<ValueType> &vectorModulus, scai::lama::Vector<ValueType> &vecDensity, scai::lama::Vector<ValueType> &vecVelocity)
{
    /* Modulus = pow(velocity,2) * Density */
    /* Velocity = sqrt( Modulus / Density )  */
    vecVelocity = vectorModulus / vecDensity;   /* = Modulus / Density */
    vecVelocity = lama::sqrt( vecVelocity );          /* = sqrt( Modulus / Density ) */
};

/*! \brief Get const reference to inverseDensity model parameter
 * 
 * If inverseDensity is dirty eg. because the density was modified, inverseDensity will be calculated from density.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensity()
{
    if (dirtyFlagInverseDensity) {
        HOST_PRINT(density.getDistributionPtr()->getCommunicatorPtr(), "Inverse density will be calculated from density\n");
        inverseDensity = 1 / density;
	dirtyFlagInverseDensity = false;
    }
    return (inverseDensity);
}

/*! \brief Get const reference to inverseDensity model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensity() const
{
    SCAI_ASSERT(dirtyFlagInverseDensity == false, "Inverse density has to be recalculated! ");
    return (inverseDensity);
}

/*! \brief Get const reference to density model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getDensity() const
{
    return (density);
}

/*! \brief Set density model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setDensity(scai::lama::Vector<ValueType> const &setDensity)
{
dirtyFlagPWaveModulus=true;
dirtyFlagSWaveModulus=true;
dirtyFlagInverseDensity = true; // If density will be changed, the inverse has to be refreshed if it is accessed
dirtyFlagAveraging = true; // If inverseDensity will be changed, averaging needs to be redone
density=setDensity;
}


/*! \brief Get const reference to first Lame model parameter
 * 
 * If P-Wave modulus is dirty eg. because the P-Wave velocity was modified, P-Wave modulus will be calculated from density and P-Wave velocity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getPWaveModulus()
{

    // If the modulus is dirty, than recalculate
    if (dirtyFlagPWaveModulus) {	   
	HOST_PRINT(velocityP.getDistributionPtr()->getCommunicatorPtr(), "P-Wave modulus will be calculated from density and P-Wave velocity\n");
        dirtyFlagPWaveModulus = false;
        calcModulusFromVelocity(velocityP,density,pWaveModulus);
    }

    return (pWaveModulus);
}

/*! \brief Get const reference to first Lame model parameter
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getPWaveModulus() const
{
    SCAI_ASSERT(dirtyFlagPWaveModulus == false, "P-Wave Modulus has to be recalculated! ");
    return (pWaveModulus);
}

/*! \brief Get const reference to second Lame Parameter sWaveModulus
 *
 * If S-Wave modulus is dirty eg. because the S-Wave Velocity was modified, S-Wave modulus will be calculated from density and S-Wave velocity.
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulus()
{

    // If the modulus is dirty, than recalculate
    if (dirtyFlagSWaveModulus) {
	HOST_PRINT(velocityS.getDistributionPtr()->getCommunicatorPtr(), "S-Wave modulus will be calculated from density and S-Wave velocity\n");    
        calcModulusFromVelocity(velocityS,density,sWaveModulus);
	dirtyFlagSWaveModulus=false;
    }

    return (sWaveModulus);
}

/*! \brief Get const reference to second Lame Parameter sWaveModulus
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulus() const
{
    SCAI_ASSERT(dirtyFlagSWaveModulus == false, "Modulus has to be recalculated! ");
    return (sWaveModulus);
}

/*! \brief Get const reference to P-wave velocity
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityP() const
{
    return (velocityP);
}

/*! \brief Set P-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setVelocityP(scai::lama::Vector<ValueType> const &setVelocityP)
{
dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
velocityP=setVelocityP;
}

/*! \brief Get const reference to S-wave velocity
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getVelocityS() const
{
    return (velocityS);
}

/*! \brief Set S-wave velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setVelocityS(scai::lama::Vector<ValueType> const &setVelocityS)
{
dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
dirtyFlagAveraging = true; // If S-Wave velocity will be changed, averaging needs to be redone
velocityS=setVelocityS;
}

/*! \brief Get const reference to tauP
 *
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauP() const
{
    return (tauP);
}

/*! \brief Set tauP velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setTauP(scai::lama::Vector<ValueType> const &setTauP)
{
dirtyFlagPWaveModulus = true; // the modulus vector is now dirty
tauP=setTauP;
}

/*! \brief Get const reference to tauS
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauS() const
{
    return (tauS);
}

/*! \brief Set tauS velocity model parameter
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setTauS(scai::lama::Vector<ValueType> const &setTauS)
{
tauS=setTauS;
dirtyFlagSWaveModulus = true; // the modulus vector is now dirty
dirtyFlagAveraging = true; // If S-Wave velocity will be changed, averaging needs to be redone
}


//! \brief Function to set elements of a single row of x-averaging density matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_DensityAverageMatrixX(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType /*NY*/, IndexType /*NZ*/)
{

    IndexType RowNumber_plusOne = rowNumber + 1;

    for (IndexType j = 0; j <= 1; j++) {
        // Insert 1.0/2.0 on diagonal elements exept for last element in every NX x NX matrix. Here: insert 1
        if (j == 0) {
            csrJALocal[countJA] = rowNumber;
            if (RowNumber_plusOne % NX == 0) {
                csrvaluesLocal[countJA] = 1.0;
            } else {
                csrvaluesLocal[countJA] = 1.0 / 2.0;
            }
            countJA++;
        } else {
            // Set elaments right to diagonal to 1.0/2.0 exept matrix elements with condition (row%NX)!=0)
            if (RowNumber_plusOne % NX != 0) {
                csrJALocal[countJA] = rowNumber + 1;
                csrvaluesLocal[countJA] = 1.0 / 2.0;
                countJA++;
            }
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Function to set elements of a single row of y-averaging density matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_DensityAverageMatrixY(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType NXNY = NX * NY;

    for (IndexType j = 0; j <= 1; j++) {
        // Insert 1.0/2.0 on diagonal elements exept for matrix elements in the last NX x NX matrix in every submatrix NXNY x NXNY. Here: insert 1
        if (j == 0) {
            csrJALocal[countJA] = rowNumber;
            if ((rowNumber % NXNY) >= (NX * (NY - 1))) {
                csrvaluesLocal[countJA] = 1.0;
            } else {
                csrvaluesLocal[countJA] = 1.0 / 2.0;
            }
            countJA++;
        } else {
            // Set elaments NX right to diagonal to 1.0/2.0 not if (row%NXNY) element of {0,...,Nx-1}
            if ((rowNumber % NXNY) <= (NX * (NY - 1) - 1)) {
                csrJALocal[countJA] = rowNumber + NX;
                csrvaluesLocal[countJA] = 1.0 / 2.0;
                countJA++;
            }
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Function to set elements of a single row of x-averaging density matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_DensityAverageMatrixZ(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ)
{

    IndexType NXNY = NX * NY;

    for (IndexType j = 0; j <= 1; j++) {
        // Insert 1.0/2.0 on diagonal elements exept for matrix elements in the last NXNY x NXNY matrix. Here: insert 1
        if (j == 0) {
            csrJALocal[countJA] = rowNumber;
            if (((rowNumber) >= (NXNY * (NZ - 1)))) {
                csrvaluesLocal[countJA] = 1.0;
            } else {
                csrvaluesLocal[countJA] = 1.0 / 2.0;
            }
            countJA++;
        } else {
            // Set elaments NXNY right to diagonal to 1.0/2.0
            if (rowNumber <= NXNY * (NZ - 1) - 1) {
                csrJALocal[countJA] = rowNumber + NXNY;
                csrvaluesLocal[countJA] = 1.0 / 2.0;
                countJA++;
            }
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Function to set elements of a single row averaging s-wave modulus matrix in xy-plane
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixXY(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType NXNY = NX * NY;

    for (IndexType j = 0; j <= 3; j++) {
        if (j <= 1) {
            if ((rowNumber % NXNY) >= (NX * (NY - 1))) {
                // Set diagonal elements in last submatrix
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j == 0) {
                        csrJALocal[countJA] = rowNumber;
                        csrvaluesLocal[countJA] = 1.0;
                        countJA++;
                    }
                } else {
                    csrJALocal[countJA] = rowNumber + j;
                    csrvaluesLocal[countJA] = 1.0 / 2.0;
                    countJA++;
                }
            } else {
                // Set diagonal elements and elements right next do diagonal elements
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j == 0) {
                        csrJALocal[countJA] = rowNumber;
                        csrvaluesLocal[countJA] = 1.0 / 2.0;
                        countJA++;
                    }
                } else {
                    csrJALocal[countJA] = rowNumber + j;
                    csrvaluesLocal[countJA] = 1.0 / 4.0;
                    countJA++;
                }
            }
        } else {
            if ((rowNumber % NXNY) <= (NX * (NY - 1) - 1)) {
                // Set elements NX and NX+1 right to diagonal elements
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j == 2) {
                        //  set last element in the submatrices
                        csrJALocal[countJA] = rowNumber + NX;
                        csrvaluesLocal[countJA] = 1.0 / 2.0;
                        countJA++;
                    }
                } else {
                    //set the other elements
                    csrJALocal[countJA] = rowNumber + NX + j - 2;
                    csrvaluesLocal[countJA] = 1.0 / 4.0;
                    countJA++;
                }
            }
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Function to set elements of a single row averaging s-wave modulus matrix in xz-plane
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixXZ(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ)
{

    IndexType NXNY = NX * NY;

    for (IndexType j = 0; j <= 3; j++) {
        if (j <= 1) {
            if (rowNumber >= (NXNY * (NZ - 1))) {
                // Set diagonal elements and elements next to diagonal in last submatrix
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j == 0) {
                        csrJALocal[countJA] = rowNumber;
                        csrvaluesLocal[countJA] = 1.0;
                        countJA++;
                    }
                } else {
                    csrJALocal[countJA] = rowNumber + j;
                    csrvaluesLocal[countJA] = 1.0 / 2.0;
                    countJA++;
                }
            } else {
                // Set diagonal elements and elements next to diagonal elements
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j == 0) {
                        csrJALocal[countJA] = rowNumber;
                        csrvaluesLocal[countJA] = 1.0 / 2.0;
                        countJA++;
                    }
                } else {
                    csrJALocal[countJA] = rowNumber + j;
                    csrvaluesLocal[countJA] = 1.0 / 4.0;
                    countJA++;
                }
            }
        } else {
            if (rowNumber <= (NXNY * (NZ - 1) - 1)) {
                // Set elements NXNY right to diagonal
                if ((rowNumber % NX) == (NX - 1)) {
                    if (j == 2) {
                        //  set last element in the submatrices
                        csrJALocal[countJA] = rowNumber + NXNY;
                        csrvaluesLocal[countJA] = 1.0 / 2.0;
                        countJA++;
                    }
                } else {
                    //set the other elements
                    csrJALocal[countJA] = rowNumber + NXNY + j - 2;
                    csrvaluesLocal[countJA] = 1.0 / 4.0;
                    countJA++;
                }
            }
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Function to set elements of a single row averaging s-wave modulus matrix in YZ-plane
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixYZ(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ)
{

    IndexType NXNY = NX * NY;
    IndexType NXNYNZ = NX * NY * NZ;

    for (IndexType j = 0; j <= 3; j++) {
        if (j <= 1) {
            if (rowNumber >= (NXNY * (NZ - 1))) {
                // Set elements in last submatrix NX x NX
                if (rowNumber >= NXNYNZ - NX) {
                    if (j == 0) {
                        csrJALocal[countJA] = rowNumber;
                        csrvaluesLocal[countJA] = 1.0;
                        countJA++;
                    }
                } else {
                    // Set elements in last submatrix NXNY x NXNY
                    csrJALocal[countJA] = rowNumber + j * NX;
                    csrvaluesLocal[countJA] = 1.0 / 2.0;
                    countJA++;
                }
            } else {
                // Set elements in diagonal submatrices NXNY x NXNY
                if ((rowNumber % NXNY) >= (NX * (NY - 1))) {
                    if (j == 0) {
                        csrJALocal[countJA] = rowNumber;
                        csrvaluesLocal[countJA] = 1.0 / 2.0;
                        countJA++;
                    }
                } else {
                    csrJALocal[countJA] = rowNumber + (j * NX);
                    csrvaluesLocal[countJA] = 1.0 / 4.0;
                    countJA++;
                }
            }
        } else {
            if (rowNumber <= (NXNY * (NZ - 1) - 1)) {
                // Set elements NXNY and NYNX+NX right to diagonal elements
                if ((rowNumber % NXNY) >= (NX * (NY - 1))) {
                    if (j == 2) {
                        //  set last element in the submatrices
                        csrJALocal[countJA] = rowNumber + NXNY;
                        csrvaluesLocal[countJA] = 1.0 / 2.0;
                        countJA++;
                    }
                } else {
                    //set the other elements
                    csrJALocal[countJA] = rowNumber + NXNY + (j - 2) * NX;
                    csrvaluesLocal[countJA] = 1.0 / 4.0;
                    countJA++;
                }
            }
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Calculate of averaging derivative matrix
/*!
 *
 \param Av Averaging Matrix Av
 \param calcNumberRowElements Member-Function to calculate number of elements in a single row
 \param setRowElements Member-Function to set the elements in a single row
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcAverageMatrix(scai::lama::Matrix<ValueType> &Av, calcNumberRowElements_AvPtr calcNumberRowElements, setRowElements_AvPtr setRowElements, IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices); //here the local indices of each process are retrieved and stored in localIndices

    /* Number of grid points */
    IndexType N = NX * NY * NZ;

    IndexType numLocalIndices = localIndices.size(); // Number of local indices
    IndexType numLocalValues = 0;                    // Number of local values of Matrix Df

    /* Calculate the number of values in each matrix */
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    IndexType read_localIndices_temp;                             // Temporary storage of the local index for the ongoing iterations

    for (IndexType i = 0; i < numLocalIndices; i++) {

        read_localIndices_temp = read_localIndices[i];

        /* Check for elements of Av */
        numLocalValues += (this->*calcNumberRowElements)(read_localIndices_temp, NX, NY, NZ);
    }

    /* Allocate local part to create local CSR storage*/
    hmemo::HArray<ValueType> valuesLocal(numLocalValues);
    hmemo::HArray<IndexType> csrJALocal(numLocalValues);
    hmemo::HArray<IndexType> csrIALocal(numLocalIndices + 1);

    /* Get WriteAccess to local part */
    hmemo::WriteAccess<IndexType> write_csrJALocal(csrJALocal);
    hmemo::WriteAccess<IndexType> write_csrIALocal(csrIALocal);
    hmemo::WriteAccess<ValueType> write_valuesLocal(valuesLocal);

    /* Set some counters to create the CSR Storage */
    IndexType countJA = 0;
    IndexType countIA = 0;
    write_csrIALocal[0] = 0;
    countIA++;

    /* Set the values into the indice arrays and the value array for CSR matrix build */
    for (IndexType i = 0; i < numLocalIndices; i++) {

        read_localIndices_temp = read_localIndices[i];

        // write AveragingMatrix
        (this->*setRowElements)(read_localIndices_temp, countJA, countIA, write_csrJALocal, write_csrIALocal, write_valuesLocal, NX, NY, NZ);
    }

    /* Release all read and write access */
    read_localIndices.release();

    write_csrJALocal.release();
    write_csrIALocal.release();
    write_valuesLocal.release();

    /* Create local CSR storage of Matrix D, than create distributed CSR matrix D */
    lama::CSRStorage<ValueType> Av_LocalCSR(numLocalIndices, N, csrIALocal, csrJALocal, valuesLocal);
    Av_LocalCSR.compress();
    Av.assignDistribute(Av_LocalCSR, dist, dist);
}

//! \brief Calculate number of row elements for density averaging matrix in x-direction
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixX(IndexType rowNumber, IndexType NX, IndexType /*NY*/, IndexType /*NZ*/)
{

    IndexType counter = 0;

    for (IndexType j = 0; j <= 1; j++) {
        if ((rowNumber % NX) == (NX - 1)) {
            if (j == 0) {
                counter++;
            }
        } else {
            counter++;
        }
    }
    return (counter);
}

//! \brief Calculate number of row elements for density averaging matrix in y-direction
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixY(IndexType rowNumber, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType counter = 0;
    IndexType NXNY = NX * NY;

    for (IndexType j = 0; j <= 1; j++) {
        if ((rowNumber % NXNY) >= (NX * (NY - 1))) {
            if (j == 0) {
                counter++;
            }
        } else {
            counter++;
        }
    }
    return (counter);
}

//! \brief Calculate number of row elements for density averaging matrix in z-direction
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixZ(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ)
{

    IndexType counter = 0;
    IndexType NXNY = NX * NY;

    for (IndexType j = 0; j <= 1; j++) {
        if (rowNumber >= (NXNY * (NZ - 1))) {
            if (j == 0) {
                counter++;
            }
        } else {
            counter++;
        }
    }
    return (counter);
}

//! \brief Calculate number of row elements for s-wave modulus averaging matrix in xy-plane
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixXY(IndexType rowNumber, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType counter = 0;
    IndexType NXNY = NX * NY;

    for (IndexType j = 0; j <= 3; j++) {
        if (j <= 1) {
            if (j == 0) {
                counter++;
            } else {
                if ((rowNumber % NX) <= ((NX - 1) - 1)) {
                    counter++;
                }
            }
        } else {
            if ((rowNumber % NXNY) <= (NX * (NY - 1) - 1)) {
                if (j == 2) {
                    counter++;
                } else {
                    if ((rowNumber % NX) != (NX - 1)) {
                        counter++;
                    }
                }
            }
        }
    }
    return (counter);
}

//! \brief Calculate number of row elements for s-wave modulus averaging matrix in xz-plane
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixXZ(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ)
{

    IndexType counter = 0;
    IndexType NXNY = NX * NY;

    for (IndexType j = 0; j <= 3; j++) {
        if (j <= 1) {
            if (j == 0) {
                counter++;
            } else {
                if ((rowNumber % NX) <= ((NX - 1) - 1)) {
                    counter++;
                }
            }
        } else {
            if (rowNumber <= (NXNY * (NZ - 1) - 1)) {
                if (j == 2) {
                    counter++;
                } else {
                    if ((rowNumber % NX) != (NX - 1)) {
                        counter++;
                    }
                }
            }
        }
    }
    return (counter);
}

//! \brief Calculate number of row elements for s-wave modulus averaging matrix in yz-plane
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template <typename ValueType>
IndexType KITGPI::Modelparameter::Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixYZ(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ)
{

    IndexType counter = 0;
    IndexType NXNY = NX * NY;

    for (IndexType j = 0; j <= 3; j++) {
        if (j <= 1) {
            if (j == 0) {
                counter++;
            } else {
                if ((rowNumber % NXNY) <= (NX * (NY - 1) - 1)) {
                    counter++;
                }
            }
        } else {
            if (rowNumber <= (NXNY * (NZ - 1) - 1)) {
                if (j == 2) {
                    counter++;
                } else {
                    if ((rowNumber % NXNY) <= (NX * (NY - 1) - 1)) {
                        counter++;
                    }
                }
            }
        }
    }
    return (counter);
}

//! \brief Calculate density averaging matrix in x-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcDensityAverageMatrixX(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcAverageMatrix(DensityAverageMatrixX, &Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixX, &Modelparameter<ValueType>::setRowElements_DensityAverageMatrixX, NX, NY, NZ, dist);
}

//! \brief Calculate density averaging matrix in y-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcDensityAverageMatrixY(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcAverageMatrix(DensityAverageMatrixY, &Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixY, &Modelparameter<ValueType>::setRowElements_DensityAverageMatrixY, NX, NY, NZ, dist);
}

//! \brief Calculate density averaging matrix in z-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcDensityAverageMatrixZ(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcAverageMatrix(DensityAverageMatrixZ, &Modelparameter<ValueType>::calcNumberRowElements_DensityAverageMatrixZ, &Modelparameter<ValueType>::setRowElements_DensityAverageMatrixZ, NX, NY, NZ, dist);
}

//! \brief Calculate s-wave modulus averaging matrix in x-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcSWaveModulusAverageMatrixXY(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcAverageMatrix(sWaveModulusAverageMatrixXY, &Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixXY, &Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixXY, NX, NY, NZ, dist);
}

//! \brief Calculate s-wave modulus averaging matrix in y-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcSWaveModulusAverageMatrixXZ(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcAverageMatrix(sWaveModulusAverageMatrixXZ, &Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixXZ, &Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixXZ, NX, NY, NZ, dist);
}

//! \brief Calculate s-wave modulus averaging matrix in z-direction
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calcSWaveModulusAverageMatrixYZ(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcAverageMatrix(sWaveModulusAverageMatrixYZ, &Modelparameter<ValueType>::calcNumberRowElements_SWaveModulusAverageMatrixYZ, &Modelparameter<ValueType>::setRowElements_SWaveModulusAverageMatrixYZ, NX, NY, NZ, dist);
}

/*! \brief calculate averaged inverse density modulus
 *
 \param vecDensity Density vector.
 \param vecInverseAvDensity Averaged inverse density vector which is calculated
 \param avDensityMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateInverseAveragedDensity(scai::lama::DenseVector<ValueType> &vecDensity, scai::lama::DenseVector<ValueType> &vecInverseAvDensity, scai::lama::Matrix<ValueType> &avDensityMatrix)
{
    Common::searchAndReplace<ValueType,IndexType>(vecDensity,0.001,0.001,1); // find vacuum grid points
    
    vecInverseAvDensity = avDensityMatrix * vecDensity;
    vecInverseAvDensity = 1 / vecInverseAvDensity;
    
    Common::searchAndReplace<ValueType,IndexType>(vecInverseAvDensity,999,0.0,2); // restore vacuum for IVF
}

/*! \brief calculate averaged s-wave modulus
 *
 \param vecSWaveModulus s-wave modulus vector
 \param vecAvSWaveModulus Averaged s-wave modulus vector which is calculated
 \param avSWaveModulusMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateAveragedSWaveModulus(scai::lama::DenseVector<ValueType> &vecSWaveModulus, scai::lama::DenseVector<ValueType> &vecAvSWaveModulus, scai::lama::Matrix<ValueType> &avSWaveModulusMatrix)
{
    Common::searchAndReplace<ValueType,IndexType>(vecSWaveModulus,1.0,1.0,1); // avoid division by zero in air/vacuum
    
    vecAvSWaveModulus = 1 / vecSWaveModulus;
    auto temp = lama::eval<lama::DenseVector<ValueType>>(avSWaveModulusMatrix * vecAvSWaveModulus);
    vecAvSWaveModulus = 1 / temp;
    
    Common::searchAndReplace<ValueType,IndexType>(vecAvSWaveModulus,4.0,0.0,1); // restore zeros for IVF
}

/*! \brief calculate averaged tauS
 *
 \param vecTauS TauS vector
 \param vecAvTauS Averaged tauS vector which is calculated
 \param avTauSMatrix Averaging matrix which is used to calculate averaged vector
 */
template <typename ValueType>
void KITGPI::Modelparameter::Modelparameter<ValueType>::calculateAveragedTauS(scai::lama::Vector<ValueType> &vecTauS, scai::lama::Vector<ValueType> &vecAvTauS, scai::lama::Matrix<ValueType> &avTauSMatrix)
{
    vecAvTauS = vecTauS;
    vecAvTauS = avTauSMatrix * vecAvTauS;
}

/*! \brief Get const reference to averaged density in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageX()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseDensityAverageX);
}

/*! \brief Get const reference to averaged density in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageY()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseDensityAverageY);
}

/*! \brief Get const reference to averaged density in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageZ()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (inverseDensityAverageZ);
}

/*! \brief Get const reference to averaged s-wave modulus in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXY()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (sWaveModulusAverageXY);
}

/*! \brief Get const reference to averaged s-wave modulus in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXZ()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (sWaveModulusAverageXZ);
}

/*! \brief Get const reference to averaged s-wave modulus in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageYZ()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (sWaveModulusAverageYZ);
}

/*! \brief Get const reference to averaged tauS in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXY()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauSAverageXY);
}

/*! \brief Get const reference to averaged tauS in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXZ()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauSAverageXZ);
}

/*! \brief Get const reference to averaged tauS in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageYZ()
{
    // If Averaging is outdated or has to be calculated for the first time, than recalculate averaging
    if (dirtyFlagAveraging == true) {
        calculateAveraging();
    }
    return (tauSAverageYZ);
}

/*! \brief Get const reference to averaged density in x-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageX() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseDensityAverageX);
}

/*! \brief Get const reference to averaged density in y-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseDensityAverageY);
}

/*! \brief Get const reference to averaged density in z-direction
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getInverseDensityAverageZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (inverseDensityAverageZ);
}

/*! \brief Get const reference to averaged s-wave modulus in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (sWaveModulusAverageXY);
}

/*! \brief Get const reference to averaged s-wave modulus in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageXZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (sWaveModulusAverageXZ);
}

/*! \brief Get const reference to averaged s-wave modulus in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getSWaveModulusAverageYZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (sWaveModulusAverageYZ);
}

/*! \brief Get const reference to averaged tauS in xy-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXY() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauSAverageXY);
}

/*! \brief Get const reference to averaged tauS in xz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageXZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauSAverageXZ);
}

/*! \brief Get const reference to averaged tauS in yz-plane
 */
template <typename ValueType>
scai::lama::Vector<ValueType> const &KITGPI::Modelparameter::Modelparameter<ValueType>::getTauSAverageYZ() const
{
    SCAI_ASSERT(dirtyFlagAveraging == false, "Averaging has to be recalculated! ");
    return (tauSAverageYZ);
}

/*! \brief Getter method for dirtyFlagPWaveModulus */
template <typename ValueType>
bool KITGPI::Modelparameter::Modelparameter<ValueType>::getDirtyFlagPWaveModulus() const
{
    return (dirtyFlagPWaveModulus);
}

/*! \brief Getter method for dirtyFlagSWaveModulus */
template <typename ValueType>
bool KITGPI::Modelparameter::Modelparameter<ValueType>::getDirtyFlagSWaveModulus() const
{
    return (dirtyFlagSWaveModulus);
}

/*! \brief Getter method for dirtyFlagInverseDensity */
template <typename ValueType>
bool KITGPI::Modelparameter::Modelparameter<ValueType>::getDirtyFlagInverseDensity() const
{
    return (dirtyFlagInverseDensity);
}

/*! \brief Getter method for dirtyFlagAveraging */
template <typename ValueType>
bool KITGPI::Modelparameter::Modelparameter<ValueType>::getDirtyFlagAveraging() const
{
    return (dirtyFlagAveraging);
}

/*! \brief Overloading = Operation
 *
 \param rhs Model which is copied.
 */
template <typename ValueType>
KITGPI::Modelparameter::Modelparameter<ValueType> &KITGPI::Modelparameter::Modelparameter<ValueType>::operator=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    assign(rhs);
    return *this;
}




/*! \brief Overloading += Operation
 *
 \param rhs Model which is added.
 */
template <typename ValueType>
KITGPI::Modelparameter::Modelparameter<ValueType> &KITGPI::Modelparameter::Modelparameter<ValueType>::operator+=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    plusAssign(rhs);
    return *this;
}

/*! \brief Overloading -= Operation
 *
 \param rhs Model which is substracted.
 */
template <typename ValueType>
KITGPI::Modelparameter::Modelparameter<ValueType> &KITGPI::Modelparameter::Modelparameter<ValueType>::operator-=(KITGPI::Modelparameter::Modelparameter<ValueType> const &rhs)
{
    minusAssign(rhs);
    return *this;
}

template class KITGPI::Modelparameter::Modelparameter<float>;
template class KITGPI::Modelparameter::Modelparameter<double>;

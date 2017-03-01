#include "Derivatives.hpp"
using namespace scai;

//! \brief Function to set elements of a single row of DybVelocity matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param read_FDCoeff_f FD-coefficients for reading
 \param read_FDCoeff_b FD-coefficients for reading
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setRowElements_DybVelocity(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType NXNY = NX * NY;
    IndexType rowEffective = rowNumber % NXNY;
    IndexType coeffPosEffective;
    IndexType coeffPosEffectiveImag;

    //Check if grid point (j/2-1) steps backward is available
    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {

        coeffPosEffective = (rowEffective - (j / 2) * NX);

        if (coeffPosEffective >= 0) {

            csrJALocal[countJA] = rowNumber - (j / 2) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_b[(j / 2 - 1)];

            /* Check for stencil outside matrix for imaging */
            coeffPosEffectiveImag = (rowEffective - ((j + 2) / 2) * NX);
            if (j + 2 <= spatialFDorder && coeffPosEffectiveImag < 0) {
                csrvaluesLocal[countJA] -= read_FDCoeff_b[(((j + 2) / 2) - 1)];
            }

            countJA++;
        }
    }

    //Check if grid point j/2 steps forward is available
    for (IndexType j = 2; j <= spatialFDorder; j += 2) {

        coeffPosEffective = rowEffective + NX * (j / 2 - 1);

        if (coeffPosEffective < NXNY) {

            csrJALocal[countJA] = rowNumber + (j / 2 - 1) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_f[j / 2 - 1];

            /* Check for stencil outside matrix for imaging */
            coeffPosEffectiveImag = (rowEffective - ((j) / 2) * NX);
            if (j <= spatialFDorder && coeffPosEffectiveImag < 0) {
                csrvaluesLocal[countJA] -= read_FDCoeff_b[(((j) / 2) - 1)];
            }

            countJA++;
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Function to set elements of a single row of DybPressure matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param read_FDCoeff_f FD-coefficients for reading
 \param read_FDCoeff_b FD-coefficients for reading
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setRowElements_DybPressure(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType NXNY = NX * NY;
    IndexType rowEffective = rowNumber % NXNY;
    IndexType coeffPosEffective;

    //Check if grid point (j/2-1) steps backward is available
    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {

        coeffPosEffective = (rowEffective - (j / 2) * NX);

        if (coeffPosEffective >= 0) {

            csrJALocal[countJA] = rowNumber - (j / 2) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_b[(j / 2 - 1)];

            /* Zeroing elements located at the feee surface */
            if (rowEffective < NX) {
                csrvaluesLocal[countJA] = 0.0;
            }

            countJA++;
        }
    }

    //Check if grid point j/2 steps forward is available
    for (IndexType j = 2; j <= spatialFDorder; j += 2) {

        coeffPosEffective = rowEffective + NX * (j / 2 - 1);

        if (coeffPosEffective < NXNY) {

            csrJALocal[countJA] = rowNumber + (j / 2 - 1) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_f[j / 2 - 1];

            /* Zeroing elements located at the feee surface */
            if (rowEffective < NX) {
                csrvaluesLocal[countJA] = 0.0;
            }

            countJA++;
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Function to set elements of a single row of Dyb matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param read_FDCoeff_f FD-coefficients for reading
 \param read_FDCoeff_b FD-coefficients for reading
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setRowElements_Dyb(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType NXNY = NX * NY;
    IndexType rowEffective = rowNumber % NXNY;
    IndexType coeffPosEffective;

    //Check if grid point (j/2-1) steps backward is available
    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {

        coeffPosEffective = (rowEffective - (j / 2) * NX);

        if (coeffPosEffective >= 0) {

            csrJALocal[countJA] = rowNumber - (j / 2) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_b[(j / 2 - 1)];

            countJA++;
        }
    }

    //Check if grid point j/2 steps forward is available
    for (IndexType j = 2; j <= spatialFDorder; j += 2) {

        coeffPosEffective = rowEffective + NX * (j / 2 - 1);

        if (coeffPosEffective < NXNY) {

            csrJALocal[countJA] = rowNumber + (j / 2 - 1) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_f[j / 2 - 1];

            countJA++;
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Function to set elements of a single row of DzfVelocity matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param read_FDCoeff_f FD-coefficients for reading
 \param read_FDCoeff_b FD-coefficients for reading
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setRowElements_DyfVelocity(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType NXNY = NX * NY;
    IndexType rowEffective = rowNumber % NXNY;
    IndexType coeffPosEffective;
    IndexType coeffPosEffectiveImag;

    //Check if grid point (j/2-1) steps backward is available
    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {

        coeffPosEffective = (rowEffective - (j / 2 - 1) * NX);

        if (coeffPosEffective >= 0) {

            csrJALocal[countJA] = rowNumber - (j / 2 - 1) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_b[(j / 2 - 1)];

            /* Check for stencil outside matrix for imaging */
            coeffPosEffectiveImag = (rowEffective - ((j + 4) / 2 - 1) * NX);
            if (j + 4 <= spatialFDorder && coeffPosEffectiveImag < 0) {
                csrvaluesLocal[countJA] -= read_FDCoeff_b[(((j + 4) / 2) - 1)];
            }

            /* Zeroing elements located at the feee surface */
            if (coeffPosEffective < NX) {
                csrvaluesLocal[countJA] = 0.0;
            }

            countJA++;
        }
    }

    //Check if grid point j/2 steps forward is available
    for (IndexType j = 2; j <= spatialFDorder; j += 2) {

        coeffPosEffective = rowEffective + NX * (j / 2);

        if (coeffPosEffective < NXNY) {

            csrJALocal[countJA] = rowNumber + (j / 2) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_f[j / 2 - 1];

            /* Check for stencil outside matrix for imaging */
            coeffPosEffectiveImag = (rowEffective - ((j + 2) / 2) * NX);
            if (j + 2 <= spatialFDorder && coeffPosEffectiveImag < 0) {
                csrvaluesLocal[countJA] -= read_FDCoeff_b[(((j + 2) / 2) - 1)];
            }

            countJA++;
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Calculate Dxf matrix
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDxf(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcDerivativeMatrix(Dxf, &Derivatives<ValueType>::calcNumberRowElements_Dxf, &Derivatives<ValueType>::setRowElements_Dxf, NX, NY, NZ, dist);
}

//! \brief Calculate Dyf matrix
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyf(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcDerivativeMatrix(Dyf, &Derivatives<ValueType>::calcNumberRowElements_Dyf, &Derivatives<ValueType>::setRowElements_Dyf, NX, NY, NZ, dist);
}

//! \brief Calculate Dzf matrix
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDzf(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcDerivativeMatrix(Dzf, &Derivatives<ValueType>::calcNumberRowElements_Dzf, &Derivatives<ValueType>::setRowElements_Dzf, NX, NY, NZ, dist);
}

//! \brief Calculate DyfPressure matrix
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyfPressure(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcDerivativeMatrix(DyfPressure, &Derivatives<ValueType>::calcNumberRowElements_Dyf, &Derivatives<ValueType>::setRowElements_Dyf, NX, NY, NZ, dist);
}

//! \brief Calculate DybVelocity matrix
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDybVelocity(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcDerivativeMatrix(DybVelocity, &Derivatives<ValueType>::calcNumberRowElements_Dyb, &Derivatives<ValueType>::setRowElements_DybVelocity, NX, NY, NZ, dist);
}

//! \brief Calculate DybPressure matrix
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDybPressure(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcDerivativeMatrix(DybPressure, &Derivatives<ValueType>::calcNumberRowElements_Dyb, &Derivatives<ValueType>::setRowElements_DybPressure, NX, NY, NZ, dist);
}

//! \brief Calculate Dyb matrix
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyb(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcDerivativeMatrix(Dyb, &Derivatives<ValueType>::calcNumberRowElements_Dyb, &Derivatives<ValueType>::setRowElements_Dyb, NX, NY, NZ, dist);
}

//! \brief Calculate DyfVelocity matrix
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyfVelocity(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcDerivativeMatrix(DyfVelocity, &Derivatives<ValueType>::calcNumberRowElements_Dyf, &Derivatives<ValueType>::setRowElements_DyfVelocity, NX, NY, NZ, dist);
}

//! \brief Calculate of derivative matrix
/*!
 *
 \param D Derivative matrix
 \param calcNumberRowElements_D Member-Function to calculate number of elements in a single row
 \param setRowElements_D Member-Function to set the elements in a single row
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDerivativeMatrix(scai::lama::Matrix &D, calcNumberRowElements_DPtr calcNumberRowElements_D, setRowElements_DPtr setRowElements_D, IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{

    /* Get local "global" indices */
    hmemo::HArray<IndexType> localIndices;
    dist->getOwnedIndexes(localIndices); //here the local indices of each process are retrieved and stored in localIndices

    /* Do some calculations to avoid it within the for loops */
    IndexType N = NX * NY * NZ;

    IndexType numLocalIndices = localIndices.size(); // Number of local indices
    IndexType numLocalValues = 0;                    // Number of local values of Matrix Df

    /* Calculate the number of values in each matrix */
    hmemo::ReadAccess<IndexType> read_localIndices(localIndices); // Get read access to localIndices
    IndexType read_localIndices_temp;                             // Temporary storage of the local index for the ongoing iterations
    for (IndexType i = 0; i < numLocalIndices; i++) {

        read_localIndices_temp = read_localIndices[i];

        /* Check for elements of D */
        numLocalValues += (this->*calcNumberRowElements_D)(read_localIndices_temp, NX, NY, NZ);
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

    /* Get ReadAccess to FD-Coefficients */
    hmemo::ReadAccess<ValueType> read_FDCoef_f(FDCoef_f);
    hmemo::ReadAccess<ValueType> read_FDCoef_b(FDCoef_b);

    /* Set the values into the indice arrays and the value array for CSR matrix build */
    for (IndexType i = 0; i < numLocalIndices; i++) {

        read_localIndices_temp = read_localIndices[i];

        /*------------*/
        /* Matrix D */
        /*------------*/
        (this->*setRowElements_D)(read_localIndices_temp, countJA, countIA, read_FDCoef_f, read_FDCoef_b, write_csrJALocal, write_csrIALocal, write_valuesLocal, NX, NY, NZ);
    }

    /* Release all read and write access */
    read_localIndices.release();
    read_FDCoef_f.release();
    read_FDCoef_b.release();

    write_csrJALocal.release();
    write_csrIALocal.release();
    write_valuesLocal.release();

    /* Create local CSR storage of Matrix D, than create distributed CSR matrix D */
    lama::CSRStorage<ValueType> D_LocalCSR(numLocalIndices, N, numLocalValues, csrIALocal, csrJALocal, valuesLocal);
    D_LocalCSR.compress();
    D.assign(D_LocalCSR, dist, dist);
}

//! \brief Function to set elements of a single row of Dzf matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param read_FDCoeff_f FD-coefficients for reading
 \param read_FDCoeff_b FD-coefficients for reading
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setRowElements_Dzf(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType NZ)
{

    IndexType rowNumber_plusOne = rowNumber + 1;
    IndexType NXNY = NX * NY;

    //Check if grid point j/2-1 steps backward is available
    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {
        if (rowNumber_plusOne > (j / 2 - 1) * NXNY) {
            csrJALocal[countJA] = rowNumber - (j / 2 - 1) * NXNY;
            csrvaluesLocal[countJA] = read_FDCoeff_b[(j / 2 - 1)];
            countJA++;
        }
    }
    //Check if grid point j/2 steps forward is available
    for (IndexType j = 2; j <= spatialFDorder; j += 2) {
        if (rowNumber_plusOne <= NXNY * (NZ - j / 2)) {
            csrJALocal[countJA] = rowNumber + NXNY * (j / 2);
            csrvaluesLocal[countJA] = read_FDCoeff_f[j / 2 - 1];
            countJA++;
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Function to set elements of a single row of Dyf matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param read_FDCoeff_f FD-coefficients for reading
 \param read_FDCoeff_b FD-coefficients for reading
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setRowElements_Dyf(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType NXNY = NX * NY;
    IndexType rowEffective = rowNumber % NXNY;
    IndexType coeffPosEffective;

    //Check if grid point (j/2-1) steps backward is available
    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {

        coeffPosEffective = (rowEffective - (j / 2 - 1) * NX);

        if (coeffPosEffective >= 0) {

            csrJALocal[countJA] = rowNumber - (j / 2 - 1) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_b[(j / 2 - 1)];

            countJA++;
        }
    }

    //Check if grid point j/2 steps forward is available
    for (IndexType j = 2; j <= spatialFDorder; j += 2) {

        coeffPosEffective = rowEffective + NX * (j / 2);

        if (coeffPosEffective < NXNY) {

            csrJALocal[countJA] = rowNumber + (j / 2) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_f[j / 2 - 1];

            countJA++;
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Function to set elements of a single row of Dxf matrix
/*!
 *
 \param rowNumber Number of current row
 \param countJA Counter for JA Elements
 \param countIA Counter for IA Elements
 \param read_FDCoeff_f FD-coefficients for reading
 \param read_FDCoeff_b FD-coefficients for reading
 \param csrJALocal Local values of JA
 \param csrIALocal Local values of IA
 \param csrvaluesLocal Local values of the matrix content
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setRowElements_Dxf(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType /*NY*/, IndexType /*NZ*/)
{

    IndexType rowNumber_plusOne = rowNumber + 1;

    //Check if grid point (j/2 - 1) steps backward is available
    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {
        if ((rowNumber_plusOne % NX >= j / 2) || (rowNumber_plusOne % NX == 0)) {
            csrJALocal[countJA] = rowNumber - (j / 2 - 1);
            csrvaluesLocal[countJA] = read_FDCoeff_b[j / 2 - 1];
            countJA++;
        }
    }
    //Check if grid point (j/2) steps forward is available
    for (IndexType j = 2; j <= spatialFDorder; j += 2) {
        if ((rowNumber_plusOne % NX <= NX - j / 2) && (rowNumber_plusOne % NX != 0)) {
            csrJALocal[countJA] = rowNumber + j / 2;
            csrvaluesLocal[countJA] = read_FDCoeff_f[j / 2 - 1];
            countJA++;
        }
    }
    csrIALocal[countIA] = countJA;
    countIA++;
}

//! \brief Calculate number of row elements for Dxf matrix
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template <typename ValueType>
IndexType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcNumberRowElements_Dxf(IndexType rowNumber, IndexType NX, IndexType /*NY*/, IndexType /*NZ*/)
{

    rowNumber = rowNumber + 1;
    IndexType counter = 0;

    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {
        //Check if grid point (j/2 - 1) steps backward is available
        if ((rowNumber % NX >= j / 2) || (rowNumber % NX == 0)) {
            counter++;
        }
        //Check if grid point (j/2) steps forward is available
        if ((rowNumber % NX <= NX - j / 2) && (rowNumber % NX != 0)) {
            counter++;
        }
    }

    return (counter);
}

//! \brief Calculate number of row elements for Dyb matrix
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Y-direction
 \return counter Number of elements in this row
 */
template <typename ValueType>
IndexType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcNumberRowElements_Dyb(IndexType rowNumber, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType counter = 0;
    IndexType NXNY = NX * NY;
    IndexType rowEffective = rowNumber % NXNY;
    IndexType coeffPosEffective;

    //Check if grid point (j/2-1) steps backward is available
    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {

        coeffPosEffective = (rowEffective - (j / 2) * NX);

        if (coeffPosEffective >= 0) {

            counter++;
        }
    }

    //Check if grid point j/2 steps forward is available
    for (IndexType j = 2; j <= spatialFDorder; j += 2) {

        coeffPosEffective = rowEffective + NX * (j / 2 - 1);

        if (coeffPosEffective < NXNY) {

            counter++;
        }
    }
    return (counter);
}

//! \brief Calculate number of row elements for Dyf matrix
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Y-direction
 \return counter Number of elements in this row
 */
template <typename ValueType>
IndexType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcNumberRowElements_Dyf(IndexType rowNumber, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    rowNumber = rowNumber + 1;
    IndexType NXNY = NX * NY;
    IndexType counter = 0;

    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {
        //Check if grid point (j/2-1) steps backward is available
        if ((rowNumber % NXNY > (j / 2 - 1) * NX) || (rowNumber % NXNY == 0)) {
            counter++;
        }
        //Check if grid point j/2 steps forward is available
        if ((rowNumber % NXNY <= NX * (NY - j / 2)) && (rowNumber % NXNY != 0)) {
            counter++;
        }
    }

    return (counter);
}

//! \brief Calculate number of row elements for Dzf matrix
/*!
 *
 \param rowNumber Row Indice
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \return counter Number of elements in this row
 */
template <typename ValueType>
IndexType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcNumberRowElements_Dzf(IndexType rowNumber, IndexType NX, IndexType NY, IndexType NZ)
{

    rowNumber = rowNumber + 1;
    IndexType NXNY = NX * NY;
    IndexType counter = 0;

    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {
        //Check if grid point j/2-1 steps backward is available
        if (rowNumber > (j / 2 - 1) * NXNY) {
            counter++;
        }
        //Check if grid point j/2 steps forward is available
        if (rowNumber <= NXNY * (NZ - j / 2)) {
            counter++;
        }
    }

    return (counter);
}

//! \brief Getter method for the spatial FD-order
template <typename ValueType>
IndexType KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getSpatialFDorder() const
{
    return (spatialFDorder);
}

//! \brief Wrapper to support configuration
/*!
 *
 \param dist Distribution of the wavefield
 \param ctx Context
 \param config Configuration
 \param comm Communicator
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::initializeMatrices(scai::dmemo::DistributionPtr dist, scai::hmemo::ContextPtr ctx, Configuration::Configuration const &config, scai::dmemo::CommunicatorPtr comm)
{
    initializeMatrices(dist, ctx, config.get<IndexType>("NX"), config.get<IndexType>("NY"), config.get<IndexType>("NZ"), config.get<ValueType>("DH"), config.get<ValueType>("DT"), config.get<IndexType>("spatialFDorder"), comm);
}

//! \brief Getter method for derivative matrix DybPressure
template <typename ValueType>
scai::lama::Matrix const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDybPressure() const
{
    if (useFreeSurface) {
        return (DybPressure);
    } else {
        return (Dyb);
    }
}
//! \brief Getter method for derivative matrix DybVelocity
template <typename ValueType>
scai::lama::Matrix const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDybVelocity() const
{
    if (useFreeSurface) {
        return (DybVelocity);
    } else {
        return (Dyb);
    }
}

//! \brief Getter method for derivative matrix DyfPressure
template <typename ValueType>
scai::lama::Matrix const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyfPressure() const
{
    if (useFreeSurface) {
        return (DyfPressure);
    } else {
        return (Dyf);
    }
}
//! \brief Getter method for derivative matrix DyfVelocity
template <typename ValueType>
scai::lama::Matrix const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyfVelocity() const
{
    if (useFreeSurface) {
        return (DyfVelocity);
    } else {
        return (Dyf);
    }
}

//! \brief Getter method for derivative matrix Dxf
template <typename ValueType>
scai::lama::Matrix const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDxf() const
{
    return (Dxf);
}

//! \brief Getter method for derivative matrix Dyf
template <typename ValueType>
scai::lama::Matrix const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyf() const
{
    return (Dyf);
}

//! \brief Getter method for derivative matrix Dzf
template <typename ValueType>
scai::lama::Matrix const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDzf() const
{
    return (Dzf);
}

//! \brief Getter method for derivative matrix Dxb
template <typename ValueType>
scai::lama::Matrix const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDxb() const
{
    return (Dxb);
}

//! \brief Getter method for derivative matrix Dyb
template <typename ValueType>
scai::lama::Matrix const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyb() const
{
    return (Dyb);
}

//! \brief Getter method for derivative matrix Dzb
template <typename ValueType>
scai::lama::Matrix const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDzb() const
{
    return (Dzb);
}

//! \brief Set FD coefficients for each order
/*!
 *
 \param spFDo Order of spatial FD-coefficient
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setFDCoef(IndexType spFDo)
{
    FDCoef_f.resize(spFDo / 2);
    FDCoef_b.resize(spFDo / 2);
    scai::hmemo::WriteAccess<ValueType> write_FDCoef_f(FDCoef_f);
    scai::hmemo::WriteAccess<ValueType> write_FDCoef_b(FDCoef_b);

    switch (spFDo) {
    case 2:
        write_FDCoef_f[0] = 1.0;
        write_FDCoef_b[0] = -1.0;
        break;
    case 4:
        write_FDCoef_f[1] = -1.0 / 24.0;
        write_FDCoef_f[0] = 9.0 / 8.0;
        write_FDCoef_b[0] = -9.0 / 8.0;
        write_FDCoef_b[1] = 1.0 / 24.0;
        break;
    case 6:
        write_FDCoef_f[2] = 3.0 / 640.0;
        write_FDCoef_f[1] = -25.0 / 384.0;
        write_FDCoef_f[0] = 75.0 / 64.0;
        write_FDCoef_b[0] = -75.0 / 64.0;
        write_FDCoef_b[1] = 25.0 / 384.0;
        write_FDCoef_b[2] = -3.0 / 640.0;
        break;
    case 8:
        write_FDCoef_f[3] = -5.0 / 7168.0;
        write_FDCoef_f[2] = 49.0 / 5120.0;
        write_FDCoef_f[1] = -245.0 / 3072.0;
        write_FDCoef_f[0] = 1225.0 / 1024.0;
        write_FDCoef_b[0] = -1225.0 / 1024.0;
        write_FDCoef_b[1] = 245.0 / 3072.0;
        write_FDCoef_b[2] = -49.0 / 5120.0;
        write_FDCoef_b[3] = 5.0 / 7168.0;
        break;
    case 10:
        write_FDCoef_f[4] = 8756999275442633.0 / 73786976294838206464.0;
        write_FDCoef_f[3] = -8142668969129685.0 / 4611686018427387904.0;
        write_FDCoef_f[2] = 567.0 / 40960.0;
        write_FDCoef_f[1] = -735.0 / 8192.0;
        write_FDCoef_f[0] = 19845.0 / 16384.0;
        write_FDCoef_b[0] = -19845.0 / 16384.0;
        write_FDCoef_b[1] = 735.0 / 8192.0;
        write_FDCoef_b[2] = -567.0 / 40960.0;
        write_FDCoef_b[3] = 8142668969129685.0 / 4611686018427387904.0;
        write_FDCoef_b[4] = -8756999275442633.0 / 73786976294838206464.0;
        break;
    case 12:
        write_FDCoef_f[5] = -6448335830095439.0 / 295147905179352825856.0;
        write_FDCoef_f[4] = 1655620175512543.0 / 4611686018427387904.0;
        write_FDCoef_f[3] = -6842103786556949.0 / 2305843009213693952.0;
        write_FDCoef_f[2] = 628618285389933.0 / 36028797018963968.0;
        write_FDCoef_f[1] = -436540475965291.0 / 4503599627370496.0;
        write_FDCoef_f[0] = 2750204998582123.0 / 2251799813685248.0;
        write_FDCoef_b[0] = -2750204998582123.0 / 2251799813685248.0;
        write_FDCoef_b[1] = 436540475965291.0 / 4503599627370496.0;
        write_FDCoef_b[2] = -628618285389933.0 / 36028797018963968.0;
        write_FDCoef_b[3] = 6842103786556949.0 / 2305843009213693952.0;
        write_FDCoef_b[4] = -1655620175512543.0 / 4611686018427387904.0;
        write_FDCoef_b[5] = 6448335830095439.0 / 295147905179352825856.0;
        break;
    default:
        COMMON_THROWEXCEPTION(" Unkown spatialFDorder value.");
        break;
    }
    write_FDCoef_f.release();
    write_FDCoef_b.release();
}

template class KITGPI::ForwardSolver::Derivatives::Derivatives<float>;
template class KITGPI::ForwardSolver::Derivatives::Derivatives<double>;
#include "Derivatives.hpp"
using namespace scai;

//! \brief Function to set elements of a single row of DybFreeSurface matrix
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
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setRowElements_DybFreeSurface(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType NXNY = NX * NY;
    IndexType rowEffective = rowNumber % NXNY;
    /* yIndex is equal to the vertical distance to the free surface */
    IndexType yIndex = IndexType(rowEffective / NX);
    IndexType coeffPosEffective;

    //loop over spatialFDorder/2 points backward of the target row
    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {

        if (rowEffective >= (j / 2) * NX) {

            /* set (normal) FD coefficients */
            csrJALocal[countJA] = rowNumber - (j / 2) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_b[(j / 2 - 1)];

            /* Vary FD coefficients at gridpoints if they have an image point */
            if (2 * yIndex >= j / 2) {
                IndexType imageIndex = 2 * yIndex - j / 2;
                if (imageIndex < spatialFDorder / 2) {
                    csrvaluesLocal[countJA] -= read_FDCoeff_b[imageIndex];
                }
            }

            countJA++;
        }
    }

    //loop over spatialFDorder/2 points forward starting at the target row
    for (IndexType j = 2; j <= spatialFDorder; j += 2) {

        coeffPosEffective = rowEffective + NX * (j / 2 - 1);

        if (coeffPosEffective < NXNY) {

            /* set (normal) FD coefficients */
            csrJALocal[countJA] = rowNumber + (j / 2 - 1) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_f[j / 2 - 1];

            /* Vary FD coefficients at gridpoints if they have an image point */
            if (2 * yIndex + j / 2 - 1 >= 0) {
                IndexType imageIndex = 2 * yIndex + j / 2 - 1;
                if (imageIndex < spatialFDorder / 2) {
                    csrvaluesLocal[countJA] -= read_FDCoeff_b[imageIndex];
                }
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

        if (rowEffective >= (j / 2) * NX) {

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

//! \brief Function to set elements of a single row of DzfFreeSurface matrix
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
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::setRowElements_DyfFreeSurface(IndexType rowNumber, IndexType &countJA, IndexType &countIA, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_f, scai::hmemo::ReadAccess<ValueType> &read_FDCoeff_b, scai::hmemo::WriteAccess<IndexType> &csrJALocal, scai::hmemo::WriteAccess<IndexType> &csrIALocal, scai::hmemo::WriteAccess<ValueType> &csrvaluesLocal, IndexType NX, IndexType NY, IndexType /*NZ*/)
{

    IndexType NXNY = NX * NY;
    IndexType rowEffective = rowNumber % NXNY;
    /* (yIndex + 1/2)  is equal to the vertical distance to the free surface */
    IndexType yIndex = IndexType(rowEffective / NX);
    IndexType coeffPosEffective;

    //loop over spatialFDorder/2 points backward starting at the target row
    for (IndexType j = spatialFDorder; j >= 2; j -= 2) {

        if (rowEffective >= (j / 2 - 1) * NX) {

            /* set (normal) FD coefficients */
            csrJALocal[countJA] = rowNumber - (j / 2 - 1) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_b[(j / 2 - 1)];

            /* Vary FD coefficients at gridpoints if they have an image point */
            if (2 * yIndex >= j / 2 - 1) {
                IndexType imageIndex = 2 * yIndex - j / 2 + 1;
                if (imageIndex < spatialFDorder / 2) {
                    csrvaluesLocal[countJA] -= read_FDCoeff_b[imageIndex];
                }
            }

            countJA++;
        }
    }

    //loop over spatialFDorder/2 points forward of the target row
    for (IndexType j = 2; j <= spatialFDorder; j += 2) {

        coeffPosEffective = rowEffective + NX * (j / 2);

        if (coeffPosEffective < NXNY) {

            /* set (normal) FD coefficients */
            csrJALocal[countJA] = rowNumber + (j / 2) * NX;
            csrvaluesLocal[countJA] = read_FDCoeff_f[j / 2 - 1];

            /* Vary FD coefficients at gridpoints if they have an image point */
            if (2 * yIndex + j / 2 >= 0) {
                IndexType imageIndex = 2 * yIndex + j / 2;
                if (imageIndex < spatialFDorder / 2) {
                    csrvaluesLocal[countJA] -= read_FDCoeff_b[imageIndex];
                }
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
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDxf(scai::dmemo::DistributionPtr dist)
{
    // Attention: keep in mind topology NZ x NY x NX

    common::Stencil1D<ValueType> stencilId(1);
    common::Stencil3D<ValueType> stencil(stencilId, stencilId, stencilFD);
    // use dist for distribution
    Dxf.define(dist, stencil);
}
//! \brief Calculate Dxf sparse matrix
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDxf(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * spatialFDorder);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        IndexType Y = coordinate.y;
        IndexType Z = coordinate.z;
        for (IndexType j = 0; j < spatialFDorder; j++) {

            IndexType X = coordinate.x;
            X += (j - spatialFDorder / 2 + 1);

            if ((X >= 0) && (X < modelCoordinates.getNX())) {
                IndexType columnIndex = modelCoordinates.coordinate2index(X, Y, Z);
                assembly.push(ownedIndex, columnIndex, stencilFD.values()[j]);
            }
        }
    }

    DxfSparse = lama::zero<SparseFormat>(dist, dist);
    DxfSparse.fillFromAssembly(assembly);
    DxfSparse *= 1 / modelCoordinates.getDH();
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
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyf(scai::dmemo::DistributionPtr dist)
{
    common::Stencil1D<ValueType> stencilId(1);
    common::Stencil3D<ValueType> stencil(stencilId, stencilFD, stencilId);
    // use dist for distribution
    Dyf.define(dist, stencil);
}

//! \brief Calculate Dyf sparse matrix
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyf(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * spatialFDorder);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        IndexType X = coordinate.x;
        IndexType Z = coordinate.z;
        for (IndexType j = 0; j < spatialFDorder; j++) {

            IndexType Y = coordinate.y;
            Y += (j - spatialFDorder / 2 + 1);

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                IndexType columnIndex = modelCoordinates.coordinate2index(X, Y, Z);
                assembly.push(ownedIndex, columnIndex, stencilFD.values()[j]);
            }
        }
    }

    DyfSparse = lama::zero<SparseFormat>(dist, dist);
    DyfSparse.fillFromAssembly(assembly);
    DyfSparse *= 1 / modelCoordinates.getDH();
}

//! \brief Calculate Dzf matrix
/*!
 *
 \param NX Number of grid points in X-direction
 aram NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDzf(scai::dmemo::DistributionPtr dist)
{
    common::Stencil1D<ValueType> stencilId(1);
    common::Stencil3D<ValueType> stencil(stencilFD, stencilId, stencilId);
    // use dist for distribution
    Dzf.define(dist, stencil);
}

//! \brief Calculate Dzf sparse matrix
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDzf(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * spatialFDorder);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        IndexType X = coordinate.x;
        IndexType Y = coordinate.y;
        for (IndexType j = 0; j < spatialFDorder; j++) {

            IndexType Z = coordinate.z;
            Z += (j - spatialFDorder / 2 + 1);

            if ((Z >= 0) && (Z < modelCoordinates.getNZ())) {
                IndexType columnIndex = modelCoordinates.coordinate2index(X, Y, Z);
                assembly.push(ownedIndex, columnIndex, stencilFD.values()[j]);
            }
        }
    }

    DzfSparse = lama::zero<SparseFormat>(dist, dist);
    DzfSparse.fillFromAssembly(assembly);
    DzfSparse *= 1 / modelCoordinates.getDH();
}
//! \brief Calculate DyfFreeSurface matrix
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyfFreeSurface(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{
    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * spatialFDorder);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        IndexType X = coordinate.x;
        IndexType Z = coordinate.z;
        for (IndexType j = 0; j < spatialFDorder; j++) {
            IndexType Y = coordinate.y;
            Y += (j - spatialFDorder / 2 + 1);

            ValueType fdCoeff = stencilFD.values()[j];

            IndexType ImageIndex = spatialFDorder - 2 - 2 * coordinate.y - j;
            if (ImageIndex >= 0)
                fdCoeff -= stencilFD.values()[ImageIndex];

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                IndexType columnIndex = modelCoordinates.coordinate2index(X, Y, Z);
                assembly.push(ownedIndex, columnIndex, fdCoeff);
            }
        }
    }

    DyfFreeSurface = lama::zero<SparseFormat>(dist, dist);
    DyfFreeSurface.fillFromAssembly(assembly);
    DyfFreeSurface *= 1 / modelCoordinates.getDH();
}

//! \brief Calculate DyfFreeSurface matrix
/*!
 *
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDybFreeSurface(Acquisition::Coordinates<ValueType> const &modelCoordinates, scai::dmemo::DistributionPtr dist)
{

    hmemo::HArray<IndexType> ownedIndexes; // all (global) points owned by this process
    dist->getOwnedIndexes(ownedIndexes);

    lama::MatrixAssembly<ValueType> assembly;
    assembly.reserve(ownedIndexes.size() * spatialFDorder);

    for (IndexType ownedIndex : hmemo::hostReadAccess(ownedIndexes)) {

        Acquisition::coordinate3D coordinate = modelCoordinates.index2coordinate(ownedIndex);

        IndexType X = coordinate.x;
        IndexType Z = coordinate.z;
        for (IndexType j = 0; j < spatialFDorder; j++) {
            IndexType Y = coordinate.y;
            Y += (j - spatialFDorder / 2);

            ValueType fdCoeff = stencilFD.values()[j];

            IndexType ImageIndex = spatialFDorder - 1 - 2 * coordinate.y - j;
            if (ImageIndex >= 0)
                fdCoeff -= stencilFD.values()[ImageIndex];

            if ((Y >= 0) && (Y < modelCoordinates.getNY())) {
                IndexType columnIndex = modelCoordinates.coordinate2index(X, Y, Z);
                assembly.push(ownedIndex, columnIndex, fdCoeff);
            }
        }
    }

    DybFreeSurface = lama::zero<SparseFormat>(dist, dist);
    DybFreeSurface.fillFromAssembly(assembly);
    DybFreeSurface *= 1 / modelCoordinates.getDH();
}

//! \brief Calculate DybFreeSurface matrix
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDybFreeSurface(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcDerivativeMatrix(DybFreeSurface, &Derivatives<ValueType>::calcNumberRowElements_Dyb, &Derivatives<ValueType>::setRowElements_DybFreeSurface, NX, NY, NZ, dist);
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

//! \brief Calculate DyfFreeSurface matrix
/*!
 *
 \param NX Number of grid points in X-direction
 \param NY Number of grid points in Y-direction
 \param NZ Number of grid points in Z-direction
 \param dist Distribution
 */
template <typename ValueType>
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDyfFreeSurface(IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
{
    calcDerivativeMatrix(DyfFreeSurface, &Derivatives<ValueType>::calcNumberRowElements_Dyf, &Derivatives<ValueType>::setRowElements_DyfFreeSurface, NX, NY, NZ, dist);
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
void KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::calcDerivativeMatrix(scai::lama::Matrix<ValueType> &D, calcNumberRowElements_DPtr calcNumberRowElements_D, setRowElements_DPtr setRowElements_D, IndexType NX, IndexType NY, IndexType NZ, scai::dmemo::DistributionPtr dist)
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
    lama::CSRStorage<ValueType> D_LocalCSR(numLocalIndices, N, csrIALocal, csrJALocal, valuesLocal);
    D_LocalCSR.compress();
    D.assignDistribute(D_LocalCSR, dist, dist);
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

        if (rowEffective >= (j / 2 - 1) * NX) {

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

        if (rowEffective >= (j / 2) * NX) {

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

//! \brief Getter method for derivative matrix DybFreeSurface
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDybFreeSurface() const
{
    return (DybFreeSurface);
}

//! \brief Getter method for derivative matrix DyfFreeSurface
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyfFreeSurface() const
{
    return (DyfFreeSurface);
}

//! \brief Getter method for derivative matrix Dxf
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDxf() const
{
    if (useSparse)
        return (DxfSparse);
    else
        return (Dxf);
}

//! \brief Getter method for derivative matrix Dyf
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyf() const
{
    if (useSparse)
        return (DyfSparse);
    else
        return (Dyf);
}

//! \brief Getter method for derivative matrix Dzf
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDzf() const
{
    if (useSparse)
        return (DzfSparse);
    else
        return (Dzf);
}

//! \brief Getter method for derivative matrix Dxb
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDxb() const
{
    if (useSparse)
        return (DxbSparse);
    else
        return (Dxb);
}

//! \brief Getter method for derivative matrix Dyb
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDyb() const
{
    if (useSparse)
        return (DybSparse);
    else
        return (Dyb);
}

//! \brief Getter method for derivative matrix Dzb
template <typename ValueType>
scai::lama::Matrix<ValueType> const &KITGPI::ForwardSolver::Derivatives::Derivatives<ValueType>::getDzb() const
{
    if (useSparse)
        return (DzbSparse);
    else
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

    const ValueType FD2[] = {-1.0, 1.0};

    const ValueType FD4[] = {1.0 / 24.0, -9.0 / 8.0, 9.0 / 8.0, -1.0 / 24.0};

    const ValueType FD6[] = {-3.0 / 640.0, 25.0 / 384.0, -75.0 / 64.0,
                             75.0 / 64.0, -25.0 / 384.0, 3.0 / 640.0};

    const ValueType FD8[] = {5.0 / 7168.0, -49.0 / 5120.0, 245.0 / 3072.0, -1225.0 / 1024.0,
                             1225.0 / 1024.0, -245.0 / 3072.0, 49.0 / 5120.0, -5.0 / 7168.0};

    const ValueType FD10[] = {-8756999275442633.0 / 73786976294838206464.0,
                              8142668969129685.0 / 4611686018427387904.0,
                              -567.0 / 40960.0,
                              735.0 / 8192.0,
                              -19845.0 / 16384.0,
                              19845.0 / 16384.0,
                              -735.0 / 8192.0,
                              567.0 / 40960.0,
                              -8142668969129685.0 / 4611686018427387904.0,
                              8756999275442633.0 / 73786976294838206464.0};

    const ValueType FD12[] = {6448335830095439.0 / 295147905179352825856.0,
                              -1655620175512543.0 / 4611686018427387904.0,
                              6842103786556949.0 / 2305843009213693952.0,
                              -628618285389933.0 / 36028797018963968.0,
                              436540475965291.0 / 4503599627370496.0,
                              -2750204998582123.0 / 2251799813685248.0,
                              2750204998582123.0 / 2251799813685248.0,
                              -436540475965291.0 / 4503599627370496.0,
                              628618285389933.0 / 36028797018963968.0,
                              -6842103786556949.0 / 2305843009213693952.0,
                              1655620175512543.0 / 4611686018427387904.0,
                              -6448335830095439.0 / 295147905179352825856.0};

    switch (spFDo) {
    case 2:
        write_FDCoef_f[0] = 1.0;
        write_FDCoef_b[0] = -1.0;
        stencilFD = common::Stencil1D<ValueType>(2, FD2);
        break;
    case 4:
        write_FDCoef_f[1] = -1.0 / 24.0;
        write_FDCoef_f[0] = 9.0 / 8.0;
        write_FDCoef_b[0] = -9.0 / 8.0;
        write_FDCoef_b[1] = 1.0 / 24.0;
        stencilFD = common::Stencil1D<ValueType>(4, FD4);
        break;
    case 6:
        write_FDCoef_f[2] = 3.0 / 640.0;
        write_FDCoef_f[1] = -25.0 / 384.0;
        write_FDCoef_f[0] = 75.0 / 64.0;
        write_FDCoef_b[0] = -75.0 / 64.0;
        write_FDCoef_b[1] = 25.0 / 384.0;
        write_FDCoef_b[2] = -3.0 / 640.0;
        stencilFD = common::Stencil1D<ValueType>(6, FD6);
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
        stencilFD = common::Stencil1D<ValueType>(8, FD8);
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
        stencilFD = common::Stencil1D<ValueType>(10, FD10);
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
        stencilFD = common::Stencil1D<ValueType>(12, FD12);
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

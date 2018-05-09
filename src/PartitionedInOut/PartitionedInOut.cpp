#include "PartitionedInOut.hpp"

using namespace scai;
using namespace KITGPI;

/*! \brief read the distributed file and redistribut the entries
 *
 \param vec Vector which is returned
 \param filename Name of the file Block (without ".%r")
 \param dist Distribution of the vector
 */
template <typename ValueType>
void KITGPI::PartitionedInOut::PartitionedInOut<ValueType>::readFromDistributedFiles(scai::lama::Vector<ValueType> &vec, std::string const &filename, scai::dmemo::DistributionPtr dist)
{
    std::string fileNameBlockIn = filename.substr(0, filename.size() - 4) + ".%r.mtx";
    vec.readFromFile(fileNameBlockIn);
    vec.redistribute(dist);
}

/*! \brief read one file with several processors and redistribut the entries
 *
 \param vec Vector which is returned
 \param filename Name of the file
 \param dist Distribution of the vector
 */
template <typename ValueType>
void KITGPI::PartitionedInOut::PartitionedInOut<ValueType>::readFromOneFile(scai::lama::Vector<ValueType> &vec, std::string const &filename, scai::dmemo::DistributionPtr dist)
{
    vec.readFromFile(filename, "BLOCK");
    vec.redistribute(dist);
}

/*! \brief redistribute vector to block distribution and write to file block
 *
 \param vec Vector which is returned
 \param filename Name of the file
 */
template <typename ValueType>
void KITGPI::PartitionedInOut::PartitionedInOut<ValueType>::writeToDistributedFiles(scai::lama::Vector<ValueType> const &vec, std::string const &filename)
{
    std::unique_ptr<lama::Vector<ValueType>> updatePtr(vec.newVector());
    lama::Vector<ValueType> &tempVector = *updatePtr;
    tempVector = vec;
    std::string distFileName = filename.substr(0, filename.size() - 4) + ".dist.mtx";
    std::string fileNameBlockOut = filename.substr(0, filename.size() - 4) + ".%r.mtx";
    lama::PartitionIO::write(tempVector.getDistribution(), distFileName);

    dmemo::DistributionPtr rowDist(new dmemo::BlockDistribution(tempVector.size(), dmemo::Communicator::getCommunicatorPtr()));
    tempVector.redistribute(rowDist);
    tempVector.writeToFile(fileNameBlockOut);
}

template class KITGPI::PartitionedInOut::PartitionedInOut<float>;
template class KITGPI::PartitionedInOut::PartitionedInOut<double>;

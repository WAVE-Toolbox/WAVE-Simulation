import scipy.io as io


def read_Seismogram(filename):
    matinfo = io.mminfo('seismogram')
    numrow = matinfo[0]
    numcol = matinfo[1]

    matread = io.mmread('seismogram').tolil()

    data = matread.toarray()
    return data, numrow, numcol

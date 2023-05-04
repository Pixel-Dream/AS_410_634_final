def readingFrames(sequence):
    """
    Creates the 6 different reading frames of our data
    within each reading frame is a list of the triplet codons required
    :param sequence:
    :return:
    """

    RF1 = []
    RF2 = []
    RF3 = []
    RF4 = []
    RF5 = []
    RF6 = []
    for i in range(0, len(sequence), 3):
        rf1codons = sequence[i:i + 3]
        RF1.append(rf1codons)
    for i in range(1, len(sequence), 3):
        rf2codons = sequence[i:i + 3]
        RF2.append(rf2codons)
    for i in range(2, len(sequence), 3):
        rf3codons = sequence[i:i + 3]
        RF3.append(rf3codons)
    # Reverse Compliment
    complement = ''
    for i in range(0, len(sequence)):
        if sequence[i] == 'A':
            complement += 'T'
        elif sequence[i] == 'T':
            complement += 'A'
        elif sequence[i] == 'G':
            complement += 'C'
        elif sequence[i] == 'C':
            complement += 'G'
    # print(complement)
    reversecomplement = complement[::-1]
    # print(reversecomplement)
    for i in range(0, len(reversecomplement), 3):
        rf4codons = reversecomplement[i:i + 3]
        RF4.append(rf4codons)
    for i in range(1, len(reversecomplement), 3):
        rf5codons = reversecomplement[i:i + 3]
        RF5.append(rf5codons)
    for i in range(2, len(reversecomplement), 3):
        rf6codons = reversecomplement[i:i + 3]
        RF6.append(rf6codons)
    readingframes = [RF1, RF2, RF3, RF4, RF5, RF6]
    return readingframes


def ORFData(readingframes):
    startcodonReadingFrames = []
    stopcodonReadingFrames = []
    for j in range(6):
        startcodonReadingFrames.append([])
        stopcodonReadingFrames.append([])
        for k in range(len(readingframes[j])):
            if readingframes[j][k] == 'ATG':
                '''if j==0:
                    startIndex=3*k
                elif j==1:
                    startIndex=(3*k)+1
                elif j==2:
                    startIndex=(3*k)+2
                elif j==3:
                    startIndex=-1*(3*k)
                elif j==4:
                    startIndex=-1*((3*k)+1)
                elif j==5:
                    startIndex=-1*((3*k)+2)'''

                startcodonReadingFrames[j].append(k)
            if readingframes[j][k] == 'TAG' or readingframes[j][k] == 'TAA' or readingframes[j][k] == 'TGA':
                '''
                if j == 0:
                    stopIndex = 3 * k
                elif j == 1:
                    stopIndex = (3 * k) + 1
                elif j == 2:
                    stopIndex = (3 * k) + 2
                elif j == 3:
                    stopIndex = -1 * (3 * k)
                elif j == 4:
                    stopIndex = -1 * ((3 * k) + 1)
                elif j == 5:
                    stopIndex = -1 * ((3 * k) + 2)'''
                stopcodonReadingFrames[j].append(k)
    return startcodonReadingFrames, stopcodonReadingFrames


def printORFs(startData, stopData, min_aa=100):
    readingFrameORFs = []
    for j in range(len(startData)):
        readingFrameORFs.append([])
        for i in range(len(startData[j])):
            if max(stopData[j]) - startData[j][i] >= min_aa:
                nearest_stop = max(stopData[j])
                for k in range(len(stopData[j])):
                    if (startData[j][i] < stopData[j][k]) and (stopData[j][k] < nearest_stop):
                        nearest_stop = stopData[j][k]
                if nearest_stop - startData[j][i] >= min_aa:
                    readingFrameORFs[j].append((startData[j][i], nearest_stop))

    return readingFrameORFs


def sequenceORFs(ORFpairs, readingframes):
    ORFs = []
    for j in range(len(ORFpairs)):
        ORFs.append([])

        for k in range(len(ORFpairs[j])):
            ORFs[j].append(readingframes[j][ORFpairs[j][k][0]:ORFpairs[j][k][1] + 1])
    return ORFs

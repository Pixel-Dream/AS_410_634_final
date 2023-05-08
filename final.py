import sys

# Input

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        # The map for all sequences
        sequences = {}
        current_header = None
        current_seq = []

        for line in file:
            # Remove leading and trailing characters
            line = line.strip()
            # Encounter a new header
            if line.startswith('>'):
                # Record the previous header and sequence data
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                    # New sequence data
                    current_seq = []
                # New header
                current_header = line[1:]
            # Sequence data
            else:
                current_seq.append(line.upper().replace(' ', ''))

        # Record the last header and sequence data
        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences


# Read ORF

def readingFrames(sequence):
    """
    Creates the 6 different reading frames of our data
    within each reading frame is a list of the triplet codons required
    :param sequence:
    :return: 123
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
    """
    Gets every occurence of the start and stop codon in the file and places it into two different lists for 6 of the
    reading frames
    Created by Ashwin Mukund
    """
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
    """
    Compares every start codon to the first instance of the stop codon. This is an inefficient function, since it has an
    unecessary for loop for the stop codon list
    An ORF is only valid for every instance start codon with the first instnace of the stop codon, since that is how
    RNA Polymerase works
    Regardless, the function examines and appends every valid ORF's locations to a tuple

    Created by Ashwin Mukund and Haowen Zhou
    """
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
    """
    Takes the tuple of the valid ORF's for each RF and then takes the original reading frame's and slices the list based
    on the locations of the tuples to get the proper ORF's

    Created by Ashwin Mukund
    """
    ORFs = []
    for j in range(len(ORFpairs)):
        ORFs.append([])

        for k in range(len(ORFpairs[j])):
            ORFs[j].append(readingframes[j][ORFpairs[j][k][0]:ORFpairs[j][k][1] + 1])
    return ORFs


# Output

def maxLenIndex(ls):
    ind = 0
    max_len = 0
    for i in range(len(ls)):
        if len(ls[i]) > max_len:
            ind = i
            max_len = len(ls[i])
    return ind


def outputORF(ORFpairs, ORFs, max_codon=15, seq_name="Untitled", max_only=True):
    rtn_text = ""
    for i in range(len(ORFs)):
        if len(ORFs[i]) > 0:
            output_index = maxLenIndex(ORFs[i])
            if max_only:
                index_array = [output_index]
            else:
                index_array = range(len(ORFs[i]))

            for ind in index_array:
                output_text = ">" + seq_name + " | " + "FRAME = " + str(i+1)
                orf_len = (1+ORFpairs[i][ind][1]-ORFpairs[i][ind][0])*3
                k = ORFpairs[i][ind][0]
                if i == 0:
                    startIndex = 3 * k + 1
                elif i == 1:
                    startIndex = (3 * k) + 2
                elif i == 2:
                    startIndex = (3 * k) + 3
                elif i == 3:
                    startIndex = -1 * (3 * k + 1)
                elif i == 4:
                    startIndex = -1 * ((3 * k) + 2)
                elif i == 5:
                    startIndex = -1 * ((3 * k) + 3)
                else:
                    startIndex = "unknown"
                output_text += " POS = " + str(startIndex) + " LEN = " + str(orf_len) + "\n"
                counter = 0
                for j in range(len(ORFs[i][ind])):
                    output_text += ORFs[i][ind][j] + " "
                    counter += 1
                    if counter >= max_codon:
                        output_text += "\n"
                        counter = 0
                rtn_text += output_text + "\n"
                print(output_text)
    return(rtn_text)


# Main function

def main():
	print("This program will find ORFs for sequences.\n When prompted, please enter a FASTA file with 1 or more sequences.")

	file_path = input("Enter a FASTA file: ")

	min_bp = int(input("Enter minimum length in bp for ORFs: "))

	seq = read_fasta(file_path)
	# Test
	#import sys
	#sys.path.append('D:/share/divergene/script/final')
	#seq = read_fasta("D:/share/divergene/script/final/sequence.fasta")
	#min_bp = 300
	
	rf = {}
	ORF_dict = {}
	ORF_pair_dict = {}
	ORF_seq_dict = {}
	min_aa = round(min_bp/3)
	
	str = ""
	
	for i in seq.keys():
	    rf[i] = readingFrames(seq[i])
	    ORF_dict[i] = ORFData(rf[i])
	    ORF_pair_dict[i] = printORFs(ORF_dict[i][0],ORF_dict[i][1],min_aa)
	    ORF_seq_dict[i] = sequenceORFs(ORF_pair_dict[i],rf[i])
	    str += outputORF(ORF_pair_dict[i],ORF_seq_dict[i],max_codon=15,seq_name=i,max_only=False)
	
	with open('output.txt', 'w') as f:
	    f.write(str)

# Call main
if __name__ == "__main__":
    main()


def maxLenIndex(ls):
    ind = 0
    max_len = 0
    for i in range(len(ls)):
        if len(ls[i]) > max_len:
            ind = i
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

#=============================================================
#
# CS 476 Project 1 by Tyler Rose
#
#=============================================================



# function for selecting the fasta sequence inputs
def fastaSelect(fastaSelection):
    
    global temp                                 # declaring temp as a global variable
    
    if fastaSelection == 1:                     # if user selects the first sequence file
        print("\nsequenceA1 Selected\n")        # print name of which sequence file is selected
        fastaFileName = ("sequenceA1.txt")      # setting file name to the file selected by the user

    elif fastaSelection == 2:                   # if second sequence is selected
        print("\nsequenceA2 Selected\n")
        fastaFileName = ("sequenceA2.txt")      
        
    elif fastaSelection == 3:                   # if third sequence is selected
        print("\nsequenceB1 Selected\n")
        fastaFileName = ("sequenceB1.txt")     
        
    elif fastaSelection == 4:                   # if fourth sequence is selected
        print("\nsequenceB2 Selected\n")
        fastaFileName = ("sequenceB2.txt")       
        
    elif fastaSelection == 5:                   # if fifth sequence is selected
        print("\nsequenceC1 Selected\n")
        fastaFileName = ("sequenceC1.txt")       
        
    elif fastaSelection == 6:                   # if sixth sequence is selected
        print("\nsequenceC2 Selected\n")
        fastaFileName = ("sequenceC2.txt")       

    else:                                       # else input is invalid
        print("\nInvalid Input. Exiting...\n")  # print error statement
        quit()                                  # exit program
    
    
    with open(fastaFileName,'r') as file:       # open file
        next(file)                              # skips first line
        temp = file.readline()                  # Read the fasta sequence and store in temp 
        temp = temp.rstrip()                    # removing trailing characters

    return temp;                                # return the sequence (temp)


# function for selecting the matrix input
def matrixSelect(matrixSelection):
    global charArray                            # declaring global variable
    global array2D                              # declaring global variable
    
    
    # BLOSUM
    if matrixSelection == 1:                   # if first matrix is selected 
        print("\nBLOSUM Selected\n")            
        matrixfileName = ("BLOSUM62.txt")       
        
    # PAM(250)
    elif matrixSelection == 2:                 # if second matrix is selected
        print("\nPAM(n) Selected\n")
        matrixfileName = ("PAM250-scores.txt")             
     
    # Hydrophobicity
    elif matrixSelection == 3:                 # if third matrix is selected
        print("\nHydrophobicity Selected\n")
        matrixfileName = ("HP.txt")  
        
    # Nucleotide
    elif matrixSelection == 4:                 # if fourth matrix is selected            
        print("\nAAnucleoPP.txt Selected!")         
        matrixfileName = ("AAnucleoPP.txt")                                    
        
    # EX2 PAM(n) matrix files
    elif matrixSelection == 5:                 # if fifth matrix is selected (more options) 
        print("\nHere are some more PAM(n) matrix file options\n")
        fileSelect = float(input("Select File Name:\n1.  PAM5scores\n2.  PAM10scores\n3.  PAM25scores\n4.  PAM50scores\n5.  PAM100scores\n6.  PAM250scores\n7.  PAM500scores\n8.  PAM1000scores\n9.  PAM2000scores\n10. PAM3000scores\n\nYour Selection: ")) 
                
        if fileSelect == 1:                    # if first matrix is selected
            print("\nPAM5scores Selected\n")            
            matrixfileName = ("PAM5scores.txt")         
           
        elif fileSelect == 2:                  # if second matrix is selected
            print("\nPAM10scores.txt Selected\n")
            matrixfileName = ("PAM10scores.txt")
            
        elif fileSelect == 3:                  # if third matrix is selected 
            print("\nPAM25scores.txt Selected\n")
            matrixfileName = ("PAM25scores.txt")
            
        elif fileSelect == 4:                  # if fourth matrix is selected
            print("\nPAM50scores.txt Selected\n")
            matrixfileName = ("PAM50scores.txt")
            
        elif fileSelect == 5:                  # if fifth matrix is selected
            print("\nPAM100scores.txt Selected\n")
            matrixfileName = ("PAM100scores.txt")
            
        elif fileSelect == 6:                  # if sixth matrix is selected
            print("\nPAM250scores.txt Selected\n")
            matrixfileName = ("PAM250scores.txt")
            
        elif fileSelect == 7:                  # if seventh matrix is selected
            print("\nPAM500scores.txt Selected\n")
            matrixfileName = ("PAM500scores.txt")
            
        elif fileSelect == 8:                  # if eighth matrix is selected
            print("\nPAM1000scores.txt Selected\n")
            matrixfileName = ("PAM1000scores.txt")
            
        elif fileSelect == 9:                  # if nineth matrix is selected
            print("\nPAM2000scores.txt Selected\n")
            matrixfileName = ("PAM2000scores.txt")
            
        elif fileSelect == 10:                 # if tenth matrix is selected
            print("\nPAM3000scores.txt Selected\n")
            matrixfileName = ("PAM3000scores.txt")
            
        else:                                  # else input is invalid                  
            print("\n\nInvalid Input. Exiting Program...\n\n")      
            quit()                                                  
        
    else:                                      # else input is invalid                
        print("\n\nInvalid Input. Exiting Program...\n\n")      
        quit()
        
        
    with open(matrixfileName, 'r') as file:                                                         # opening selecteed file
                
                next(file)                                                                          # skip first line
                chars = file.readline()                                                             # read sequence line to get characters
                charArray = chars.split(',')                                                        # split at each comma
                charArray[len(charArray) - 1] = charArray[len(charArray) - 1].rstrip()              # store each char and remove trailing characters
                
                row = len(charArray)                                                                # get row length
                column = row                                                                        # set culumn length = to row length (since its a matrix)
                array2D = column * [[0] * row]                                                      # initializing 2-D array to store the values

                for i in range(row):                                                                # for loop to complete each row
                    data = file.readline()                                                          # read each line
                    data = data.split(',')                                                          # split at each comma
                    array2D[i] = data                                                               # fill array with data
                    array2D[i][len(array2D) - 1] = array2D[i][len(array2D) - 1].rstrip()            # store data in array and remove trailing characters

    return charArray, array2D;                                                                      # return charArray and the 2D array
    

# function for calculating Global OPT alignment
def globalFunction(gapPenalty, ONEfastaSeq, TWOfastaSeq, charArray, array2D):

    m = 1 + len(ONEfastaSeq)                                        # getting length of the first sequence
    n = 1 + len(TWOfastaSeq)                                        # getting length of the second sequence
    matrix = [[0 for i in range(m)] for j in range(n)]              # initializing a matrix for the values
    matrix[0][0] = 0                                                # initialize matrix elements to be zero
    dirCount = [['N' for i in range(m)] for j in range(n)]          # setting up the directional counter for retracing our steps, 'N' being none

    
    for i in range(1, m):                                           # for loop for filling in matrix 
        matrix[0][i] = gapPenalty + matrix[0][i - 1]                # filling in matrix for first sequence
    
    for j in range(1, n):                                           # for loop for filling in matrix
        matrix[j][0] = gapPenalty + matrix[j - 1][0]                # filling in matrix for for second sequence
        dirCount[j][0] = 'N'                                        # setting directional count to 'N' for none

        posB = TWOfastaSeq[j - 1]                                   # finding the character's position in the sequence
        idxB = charArray.index(posB)                                # indexing for the character's position

        for i in range(1, m):                                       # for loop to parse through first sequence length
            posA = ONEfastaSeq[i - 1]                               # finding the character's position in the sequence
            idxA = charArray.index(posA)                            # indexing for the character's position
            tempMatrix = int(float(array2D[idxA][idxB]))            # storing positions in a temp array

            maxData = max((gapPenalty + matrix[j - 1][i]),          # comparing directional values (vertical)
                          (gapPenalty + matrix[j][i - 1]),          # horizontal
                          (tempMatrix + matrix[j - 1][i - 1]))      # diagonal    
            matrix[j][i] = maxData                                  # store max directional values

            if maxData == gapPenalty + matrix[j - 1][i]:            # if Vertical was highest value
                dirCount[j][i] = 'V'                                # set directional count to vertical
            elif maxData == gapPenalty + matrix[j][i - 1]:          # if Horizontal was highest value
                dirCount[j][i] = 'H'                                # set directional count to horizontal
            else:                                                   # else Diagonal is highest value
                dirCount[j][i] = 'D'                                # set directional count to diagonal
    
    a = ''                                                          # initialize output
    b = ''                                                          # initialize output
    x = m - 1                                                       # initializing place holder for alignment
    y = n - 1                                                       # initializing place holder for alignment
    currentData = dirCount[y][x]                                    # initializing place holder for directional value

    while currentData != 'N':                                       # while loop for when a direction has been assigned 
        if currentData == 'V':                                      # if vertical
            a = '-' + a                                             # append - 
            b = TWOfastaSeq[y-1] + b                                # move to next char
            y = y - 1                                               # remove previous char
            
        elif currentData == 'H':                                    # if horizontal
            a = ONEfastaSeq[x-1] + a                                # move to previous char
            b = '-' + b                                             # append -
            x = x - 1                                               # remove char
            
        else:                                                       # else diagonal
            a = ONEfastaSeq[x-1] + a                                # move to previous char
            b = TWOfastaSeq[y-1] + b                                # move to previous char
            x = x - 1                                               # remove char
            y = y - 1                                               # remove char

        currentData = dirCount[y][x]                                # store sequence char
        score = matrix[n-1][m-1]                                    # store overall score value

    for i in range(0, n):                                           # print matrix
        print( matrix[i] )

    print("\nTop  OPT Alignment: " + a)                             # print top OPT alignment
    print("Side OPT Alignment: " + b)                               # print side OPT alignment
    print("OPT Score: " + str(score))                               # print score
    
    return


# function for calculating Semi-global OPT alignment
def semiGlobalFunction(gapPenalty, ONEfastaSeq, TWOfastaSeq, charArray, array2D):

    m = 1 + len(ONEfastaSeq)                                        # getting length of the first sequence
    n = 1 + len(TWOfastaSeq)                                        # getting length of the second sequence
    matrix = [[0 for i in range(m)] for j in range(n)]              # initializing a matrix for the values
    matrix[0][0] = 0                                                # initialize matrix elements to be zero
    dirCount = [['N' for i in range(m)] for j in range(n)]          # setting up the directional counter for retracing our steps, 'N' being none       

    pointA = 0                                                      # position of max score
    pointB = 0                                                      # position of max score
    scoreMax = 0                                                    # max score

    for i in range(1, m):                                           # for loop for filling in matrix
        matrix[0][i] = 0                                            # initializing matrix
        dirCount[0][i] = 'N'                                        # initialzing directional counter
            
    for j in range(1, n):                                           # for loop for filling in matrix
        matrix[j][0] = 0                                            # initializing matrix
        dirCount[j][0] = 'N'                                        # initialzing directional counter

        posB = TWOfastaSeq[j - 1]                                   # finding the character's position in the sequence
        idxB = charArray.index(posB)                                # indexing for the character's position

        for i in range(1, m):                                       # for loop to parse through first sequence length
            posA = ONEfastaSeq[i - 1]                               # finding the character's position in the sequence
            idxA = charArray.index(posA)                            # indexing for the character's position
            tempMatrix = int(float(array2D[idxA][idxB]))            # storing positions in a temp array

            maxData = max((gapPenalty + matrix[j - 1][i]),          # comparing directional values (vertical)
                          (gapPenalty + matrix[j][i - 1]),          # horizontal
                          (tempMatrix + matrix[j - 1][i - 1]))      # diagonal    
            matrix[j][i] = maxData                                  # store max directional values

            if maxData == gapPenalty + matrix[j - 1][i]:            # if Vertical was highest value
                dirCount[j][i] = 'V'                                # set directional count to vertical
            elif maxData == gapPenalty + matrix[j][i - 1]:          # if Horizontal was highest value
                dirCount[j][i] = 'H'                                # set directional count to horizontal
            else:                                                   # else Diagonal is highest value
                dirCount[j][i] = 'D'                                # set directional count to diagonal

            if matrix[j][i] > scoreMax and i == m - 1 or j == n - 1:# if max is on the edge 
                scoreMax = matrix[j][i]                             # set max score
                pointA = j                                          # get location
                pointB = i                                          # get location
                
   

    a = ''                                                          # initialize output
    b = ''                                                          # initialize output
    y = pointA                                                      # initializing place holder for alignment
    x = pointB                                                      # initializing place holder for alignment
    currentData = dirCount[y][x]                                    # initializing place holder for directional value
   
    while currentData != 'N':                                       # while loop for when a direction has been assigned 
        if currentData == 'V':                                      # if vertical
            a = '-' + a                                             # append - 
            b = TWOfastaSeq[y-1] + b                                # move to next char
            y = y - 1                                               # remove previous char
            
        elif currentData == 'H':                                    # if horizontal
            a = ONEfastaSeq[x-1] + a                                # move to previous char
            b = '-' + b                                             # append -
            x = x - 1                                               # remove char
            
        else:                                                       # else diagonal
            a = ONEfastaSeq[x-1] + a                                # move to previous char
            b = TWOfastaSeq[y-1] + b                                # move to previous char
            x = x - 1                                               # remove char
            y = y - 1                                               # remove char

        currentData = dirCount[y][x]                                # store sequence char

    for i in range(0, n):                                           # print matrix
        print( matrix[i] )

    print("\nTop  OPT Alignment: " + a)                             # print top OPT alignment
    print("Side OPT Alignment: " + b)                               # print side OPT alignment
    print("OPT Score: " + str(scoreMax))                            # print score
    
    return


# function for calculating Local OPT alignment
def localFunction(gapPenalty, ONEfastaSeq, TWOfastaSeq, charArray, array2D):
    
    m = 1 + len(ONEfastaSeq)                                        # getting length of the first sequence
    n = 1 + len(TWOfastaSeq)                                        # getting length of the second sequence
    matrix = [[0 for i in range(m)] for j in range(n)]              # initializing a matrix for the values
    matrix[0][0] = 0                                                # initialize matrix elements to be zero
    dirCount = [['N' for i in range(m)] for j in range(n)]          # setting up the directional counter for retracing our steps, 'N' being none       

    pointA = 0                                                      # position of max score
    pointB = 0                                                      # position of max score
    scoreMax = 0                                                    # max score
 
    for i in range(1, m):                                           # for loop for filling in matrix
        matrix[0][i] = 0                                            # initializing matrix
        dirCount[0][i] = 'N'                                        # initialzing directional counter
       
    for j in range( 1, n ):                                         # for loop for filling in matrix
        matrix[j][0] = 0                                            # initializing matrix
        dirCount[j][0] = 'N'                                        # initialzing directional counter

        posB = TWOfastaSeq[j - 1]                                   # finding the character's position in the sequence
        idxB = charArray.index(posB)                                # indexing for the character's position

        for i in range(1, m):                                       # for loop to parse through first sequence length
            posA = ONEfastaSeq[i - 1]                               # finding the character's position in the sequence
            idxA = charArray.index(posA)                            # indexing for the character's position
            tempMatrix = int(float(array2D[idxA][idxB]))            # storing positions in a temp array

            maxData = max(0,                                        # comparing directional values
                          (gapPenalty + matrix[j - 1][i]),          # vertical
                          (gapPenalty + matrix[j][i - 1]),          # horizontal
                          (tempMatrix + matrix[j - 1][i - 1]))      # diagonal    
            matrix[j][i] = maxData                                  # store max directional values

            if maxData == 0:
                dirCount[j][i] = 'N'                                # checking for directional value
            elif maxData == gapPenalty + matrix[j - 1][i]:          # if Vertical was highest value
                dirCount[j][i] = 'V'                                # set directional count to vertical
            elif maxData == gapPenalty + matrix[j][i - 1]:          # if Horizontal was highest value
                dirCount[j][i] = 'H'                                # set directional count to horizontal
            else:                                                   # else Diagonal is highest value
                dirCount[j][i] = 'D'                                # set directional count to diagonal
               
            if matrix[j][i] > scoreMax:                             # if data is greater than max score
                scoreMax = matrix[j][i]                             # set new max score
                pointA = j                                          # get location
                pointB = i                                          # get location

    a = ''                                                          # initialize output
    b = ''                                                          # initialize output
    y = pointA                                                      # initializing place holder for alignment
    x = pointB                                                      # initializing place holder for alignment
    currentData = dirCount[y][x]                                    # initializing place holder for directional value

    while currentData != 'N':                                       # while loop for when a direction has been assigned 
        if currentData == 'V':                                      # if vertical
            a = '-' + a                                             # append - 
            b = TWOfastaSeq[y-1] + b                                # move to next char
            y = y - 1                                               # remove previous char
            
        elif currentData == 'H':                                    # if horizontal
            a = ONEfastaSeq[x-1] + a                                # move to previous char
            b = '-' + b                                             # append -
            x = x - 1                                               # remove char
            
        else:                                                       # else diagonal
            a = ONEfastaSeq[x-1] + a                                # move to previous char
            b = TWOfastaSeq[y-1] + b                                # move to previous char
            x = x - 1                                               # remove char
            y = y - 1                                               # remove char

        currentData = dirCount[y][x]                                # store sequence char

    for i in range(0, n):                                           # print matrix
        print( matrix[i] )

    print("\nTop  OPT Alignment: " + a)                             # print top OPT alignment
    print("Side OPT Alignment: " + b)                               # print side OPT alignment
    print("OPT Score: " + str(scoreMax))                            # print score
    
    return

    
def main():

    # SELECT SEQUENCE TYPE
    typeSeqSel = float(input("\nSelect A Sequence Type:\n1.  Nucleotide\n2.  Peptide\n\nYour Selection: "))  # ask for user input

    if typeSeqSel == 1:                                             # if user selects Nucleotide (option 1)
        print("\nNucleotide Selected!")                             # print the selection
        
    elif typeSeqSel == 2:                                           # if user selects Peptide Sequence (protein or amino acid)
        print("\nPeptide Selected!")                                

    else:                                                           # else exit program
        print("\n\nInvalid Input. Exiting Program...\n\n")
        quit() 


    # SELECT FIRST FASTA INPUT FILE BY CALLING FUNCTION, WHICH RETURNS THE FASTA FILE NAME
    fastaSelectionONE = float(input("\nSelect A FASTA Filename:\n1.  sequenceA1\n2.  sequenceA2\n3.  sequenceB1\n4.  sequenceB2\n5.  sequenceC1\n6.  sequenceC2\n\nYour Selection: "))
    fastaSelect(fastaSelectionONE) 

    ONEfastaSeq = temp                                              # set temp to ONEfastaSeq
    print(ONEfastaSeq)                                              # print sequence from array


    # SELECT SECOND FASTA INPUT FILE
    fastaSelectionTWO = float(input("\nSelect A FASTA Filename:\n1.  sequenceA1\n2.  sequenceA2\n3.  sequenceB1\n4.  sequenceB2\n5.  sequenceC1\n6.  sequenceC2\n\nYour Selection: "))

    if fastaSelectionONE == fastaSelectionTWO:                      # if the same sequence is selected twice 
        print("\n\nPlease select a different sequence than the first selection. Exiting...\n\n")
        quit()                                                      # quit

    fastaSelect(fastaSelectionTWO)                                  # call function again 
    TWOfastaSeq = temp                                              # store sequence in temp
    print(TWOfastaSeq)                                              # print sequence

 
    # SELECT MATRIX TYPE
    matrixSelection = float(input("\nSelect a matrix:\n1.  BLOSUM\n2.  PAM(250)\n3.  Hydrophobicity\n4.  nAAnucleoPP\n5.  View more PAM(n) matrix file options\n\nYour Selection: ")) # recieving user input for matrixSelection
    matrixSelect(matrixSelection)                                   # call matrix function
    
    
    # SELECT ALIGNMENT TYPE
    alignSelection = float(input("\nSelect an alignment:\n1.  Global\n2.  Semi-Global\n3.  Local\n\nYour Selection: ")) # recieving user input
    
    if alignSelection == 1:                                         # if GLOBAL is selected
        print("\nGlobal Selected\n")
        gapPenalty = float(input("\nInput Desired Gap Penalty: "))
        print("\n")
        
        gapPenalty = abs(gapPenalty)                                # make input positive, to ensure we can make the gap penalty negative
        gapPenalty = gapPenalty * -1                                # make input negative
        
        globalFunction(gapPenalty, ONEfastaSeq, TWOfastaSeq, charArray, array2D)             
      
    elif alignSelection == 2:                                       # if SEM_GLOBAL is selected
        print("\nSemi-Global Selected\n")
        gapPenalty = float(input("\nInput Desired Gap Penalty: "))
        print("\n")    
        
        gapPenalty = abs(gapPenalty)   
        gapPenalty = gapPenalty * -1   
        
        semiGlobalFunction(gapPenalty, ONEfastaSeq, TWOfastaSeq, charArray, array2D)                 
        
    elif alignSelection == 3:                                       # if LOCAL is selected
        print("\nLocal Selected\n")
        gapPenalty = float(input("\nInput Desired Gap Penalty:"))
        print("\n")
        
        gapPenalty = abs(gapPenalty)   
        gapPenalty = gapPenalty * -1    
        
        localFunction(gapPenalty, ONEfastaSeq, TWOfastaSeq, charArray, array2D)
        
    else:
        print("\n\nInvalid Input. Exiting Program...\n\n")
        quit()

    print("\nAlignment Complete!\n")                                # message to let user know the program is complete 
    
    
if __name__ == '__main__':                                          # telling the program to start with the 'main' function
    main()

##############################################################
#
# CS 476 Project 2 by Tyler Rose
#
##############################################################


# UPGMA Algorithm
def upgma(nodes, letters, matrixOG):

    format = letters                            # declaring format equal to letters list
    tmpMatrix = matrixOG                        # declaring and storing the original distance matrix in a temp matrix
    tmpLetterList = letters                     # declaring a temp list for holding the letters
    cluster = nodes                             # declaring cluster and equalling it to nodes

    while(cluster > 1):                         # while loop that continues until there are no more clusters to merge

        minDistance = 1000                      # initializing minDistance to be extremely high to ensure it is greater than the max value in the matrix

        for i in range(cluster-1):              # for loop to get the number of clusters - 1
            for j in range(i+1, cluster):       # for loop to get the number of characters and clusters
                temp = float(tmpMatrix[i][j])   # declare temp variable to hold numer of values
                if temp <= minDistance:         # if temp is less than or equal to the minDistance
                    minDistance = temp          # set minDistance equal to temp
                    minI = i                    # store i location 
                    minJ = j                    # store j location  

        cluster -= 1                                                        # decrement cluster count
        phyTree = tmpLetterList.copy()                                      # copy the letter list over to phyTree 
        phyTree[minI] = phyTree[minI] + phyTree[minJ]                       # store the new merged pair in the phyTree
        phyTree.remove(phyTree[minJ])                                       # remove the old element
        minDistance = "{:.2f}".format(minDistance)                          # initialize minDistance to only have 2 decimal places
        newMatrix = [[0 for i in range(cluster)] for j in range(cluster)]   # initialize new matrix with updated number of clusters

        if cluster != 1:
            for i in range(cluster):                                        # if there are more than 1 cluster left
                for j in range(cluster):                                    # for loop for the number of clusters 
                    if i != minI and j != minI:                             # if the i and j do not match the location of the min distance
                        tmpIndex1 = tmpLetterList.index(phyTree[i])         # index to find the location
                        tmpIndex2 = tmpLetterList.index(phyTree[j])         # index to find the location
                        newMatrix[i][j] = tmpMatrix[tmpIndex1][tmpIndex2]   # update the new matrix with the correct location
                    elif i == j:                                            # else if i and j are equal, just keep chugging along
                        continue                                            # why are you stopping? nothing to see here. keep it moving
                    elif i == minI or j == minI:                            # else if one of the locations is correct
                        if i == minI:                                       # if i equalled minI
                            letter = phyTree[j]                             # set letter equal to j character
                        else:                                               
                            letter = phyTree[i]                             # else set letter equal to i character

                        value1 = tmpMatrix[minI][tmpLetterList.index(letter)] * len(tmpLetterList[minI])    # store value for calculations
                        value2 = tmpMatrix[minJ][tmpLetterList.index(letter)] * len(tmpLetterList[minJ])    # store value for calculation
                        newMatrix[i][j] = (float(value1) + float(value2)) / (len(tmpLetterList[minI]) + len(tmpLetterList[minJ]))   # do calculations and store answer in new matrix 

        print(f"{'Merging ' + tmpLetterList[minI] + ' & ' + tmpLetterList[minJ] + ':' : <20}", end = '   ')     # printing the clusters being merged
        print(f"{'Distance = ' + minDistance : <18}", end = '  --->  ')                                         # printing the distance
        print(phyTree, "\n")                                                                                    # print updated phylogenetic tree

        format[minI] = '(' + format[minI] + ',' + format[minJ] + ')'        # formatting the phylogenetic tree for Newick format
        format.remove(format[minJ])                                         # formatting still
        
        tmpLetterList = phyTree                                             # reset the temp letter list to the updated tree
        tmpMatrix = newMatrix                                               # reset the temp matrix to the updated matrix

    print('\nValid Phylogenetic Tree: ' + format[0])                        # print final phylogenetic tree in Newick format


# Neighbor Join Algorithm
def neighborJoin(nodes, letters, matrixOG):
    
    r = [0 for i in range(nodes)]                                           # initializing r values
    format = letters                                                        # declaring format equal to letters
    tmpMatrix = matrixOG                                                    # declaring a temp matrix equal to the original matrix
    tmpLetterList = letters.copy()                                          # copying over the letters list to a temp list
    cluster = nodes                                                         # initialzing clusters equal to nodes

    while(cluster > 2):                                                     # while loop ensuring clusters are greater than 2
        for i in range(cluster):                                            # for loop for the number of clusters
            temp = 0                                                        # initializing temp int
            for j in range(cluster):                                        # for loop counting clusters
                temp = temp + int(tmpMatrix[i][j])                          # set temp equal to the value of the matrix
            r[i] = temp / (cluster - 2)                                     # calculate the r value and store it in r array


        print('\nCalculated R Values:\n')                                       # print the r values
        for i in range(cluster):                                                # for loop for printing the r values
            print('    ' + tmpLetterList[i] + '=> ' + "{:.2f}".format(r[i]))    # print r values to 2 decimal places

        minDistance = 1000                                                      # initialize min distance to ensure its greater than the largest distance in the matrix

        print('\n\nTransition Distance Matrix:\n')                              # print transition matrix
        
        for i in range(cluster-1):                                              # for loop counting clusters
            for j in range (i+1, cluster):                                      # for loop getting proper range
                tranDist = int(tmpMatrix[i][j]) - r[i] - r[j]                   # equation for the transition distances
                print( '   TD of (' + tmpLetterList[i] + ', ' + tmpLetterList[j] + ') = ' + "{:.2f}".format(tranDist))  # print distances
                if tranDist < minDistance:                                      # if statement to find min transition distance value
                    minDistance = tranDist                                      # set new min value
                    minI = i                                                    # set i location
                    minJ = j                                                    # set j location

        br1 = (int(tmpMatrix[i][j]) + r[minI] - r[minJ]) / 2                    # calculation for branch distances
        br1 = "{:.2f}".format(br1)                                              # formatting the branch distance to 2 decimal places
        br2 = (int(tmpMatrix[i][j]) + r[minJ] - r[minI]) / 2                    # calculating for branch distances
        br2 = "{:.2f}".format(br2)                                              # formatting the branch distance to 2 decimal places


        print('\n\nClusters Merging: ' + tmpLetterList[minI] + ' & ' + tmpLetterList[minJ] + '\n')  # print the clusters that are merging
        print('   ' + tmpLetterList[minI] + ' distance = ' + br1)                                   # printing branch distance
        print('   ' + tmpLetterList[minJ] + ' distance = ' + br2 + '\n')                            # printing branch distance

        for i in range(cluster):                                                                    # for loop for cluster number
            if tmpLetterList[minI] != tmpLetterList[i] and tmpLetterList[minJ] != tmpLetterList[i]: # if current letters are not merging
                distance1 = int(tmpMatrix[i][minI])                                                 # find min distances
                distance2 = int(tmpMatrix[i][minJ])                                                 # find min distances

                tmpMatrix[minI][i] = (distance1 + distance2 - int(tmpMatrix[minI][minJ])) / 2       # calculate new distance
                tmpMatrix[i][minI] = tmpMatrix[minI][i]                                             # store new distance in matrix


        tmpLetterList[minI] = tmpLetterList[minI] + tmpLetterList[minJ]                             # merge clusters in matrix
        tmpLetterList.remove( tmpLetterList[minJ] )                                                 # remove old cluster
        format[minI] = '(' + format[minI] + ','+format[minJ] + ')'                                  # formatting
        format.remove(format[minJ])                                                                 # remove old j values
        cluster -= 1                                                                                # decrement cluster count
        r.remove(r[minJ])                                                                           # remove old r values

    format[0] = '(' + format[0] + ',' + format[1] + ')'                                             # formatting the Newick Tree
    print('\n\nFinal Merged Cluster Distance --> ' + tmpLetterList[0] + ' & ' + tmpLetterList[1] + ' = ' + tmpMatrix[0][1] + '' + tmpMatrix[1][0] )   # printing final distance                                                  # print 
    print('\nValid Phylogenetic Tree: : ' + format[0] + '\n\n')                                     # print valid phylogenetic tree


# main
def main():

    # Request user to select an Implementation
    typeSelect = float(input("\nSelect An Implementation Type:\n\n1.  UPGMA\n2.  Neighbor Join\n\nYour Selection: "))  # ask for user input

    if typeSelect == 1:                                             # if user selects UPGMA (option 1)
        print("\nUPGMA Selected!\n")                                # print the selection
    elif typeSelect == 2:                                           # if user selects Neighbor join
        print("\nNeighbor Join Selected!\n")                        # print the selection      
    else:                                                           # else exit program
        print("\n\nInvalid Input. Exiting Program...\n\n")
        quit()
        
        
    # Request user to select an input file
    fileSelect = float(input("\nSelect A File For Input:\n1.  DM-p127.txt\n2.  DM-p139.txt\n\nYour Selection: "))  # ask for user input

    if fileSelect == 1:                                             # if user selects DM-p127.txt (option 1)
        print("\nDM-p127.txt Selected!\n\n")                        # print the selection
        fileName = ("DM-p127.txt")                                  # setting file name
        
    elif fileSelect == 2:                                           # if user selects DM-p139.txt
        print("\nDM-p139.txt Selected!\n\n")                        # print the selection    
        fileName = ("DM-p139.txt")                                  # setting file name

    else:                                                           # else exit program
        print("\n\nInvalid Input. Exiting Program...\n\n")
        quit()    

    # reading from input file
    with open(fileName, 'r') as file:                               # opening matrix file
        nodes = file.readline().strip()                             # storing nodes
        nodes = int(nodes)                                          # converting nodes to integers
        letters = file.readline()                                   # storing letters
        letters = letters.split()                                   # splitting the letters into separate elements in the list

        matrixOG = [[0 for i in range(nodes)] for j in range(nodes)]# storing original (gangsta) matrix in array

        for i in range(nodes):                                      # for loop to count all nodes             
            values = file.readline()                                # storing i values in array  
            values = values.split()                                 # splitting each value to be a separate element
            matrixOG[i] = values                                    # store i values in matrix array
    

    # printing the matrix array to the console
    printChars = ''                                                 # declare variable
    printMatrix = ''                                                # clearing the data from the matrix array
    
    for i in range(nodes):                                          # for loop to count every node
        if i == 0:                                                  # if first character
            printChars = printChars + '  '                          # add a space at beginning of array for readability
        printChars = printChars + letters[i] + ' '                  # storing the characters in the matrix array
        
    for i in range(nodes):                                          # for loop to count every node
        if i == 0:                                                  # if in first position
            printMatrix = printMatrix + letters[i] + ' '            # include characters on side of matrix
        if i != 0:                                                  # if i is not zero, print a new line
            printMatrix = printMatrix + '\n' + letters[i] + ' '     # new line 
        for j in range(nodes):                                      # for loop to count every node
            printMatrix = printMatrix + matrixOG[i][j] + ' '        # store the full matrix
    
    print("Initial Distance Matrix\n")                              # print statement
    print(printChars)                                               # print characters
    print(printMatrix)                                              # print matrix
    print("\n")                                                     # print new line for readability

    
    if typeSelect == 1:                                             # if the user selects UPGMA
        upgma(nodes, letters, matrixOG)                             # call upgma function
        
    
    elif typeSelect == 2:                                           # if the user selects neighbor join
        neighborJoin(nodes, letters, matrixOG)                      # call neighborJoin function
        
    else:                                                           # else exit program
        print("\n\nInvalid Input. Exiting Program...\n\n")  
        quit()

if __name__ == '__main__':                                          # start program at main()
    main()

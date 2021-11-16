"""
Protein Sequencing Project
Name: Rachana
Roll Number:2021501004
"""

import hw6_protein_tests as test
import numpy
import matplotlib

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    f=open(filename,"r")
    text=f.read()
    #print(text)
    f.close()
    text=text.replace("\n","")
    return text


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    #print(dna)
    dna=dna.replace("T","U")
    new_dna=[]
    result=[]
    stop_Code=["UAA","UAG","UGA"]
    for i in range(startIndex,len(dna),3):
        new_dna.append(dna[i:i+3])
    #print(new_dna)
    #if new_dna[startIndex]=="ATG":
    for each in new_dna:
        if each not in stop_Code:
             result.append(each)
        else:
            result.append(each)
            break
    
    return result

'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    f=open(filename,"r")
    result_Dict={}
    j=json.load(f)
    #j=j.replace("T","U")
    #print(j)
    for key,value in j.items():
        for each in value:
            each=each.replace("T","U")
            result_Dict[each]=key
    f.close()
    return result_Dict


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    #print(codonD)
    protein =[]
    for each in codons:
        if  each=="AUG" and "Start" not in protein:
                protein.append("Start")
        else: 
            protein.append(codonD[each]) 
    #print(protein)

    return protein 


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna_File=readFile(dnaFilename)
    codon_Dict=makeCodonDictionary(codonFilename)
    #print(dna_File,"\n",codon_Dict)
    count=0
    i=0
    final_protein_List=[]
    while i< len(dna_File):
        if dna_File[i:i+3]=="ATG":
            startIndex=i
            make_DnaToRna=dnaToRna(dna_File,startIndex)
            make_Protein=generateProtein(make_DnaToRna,codon_Dict)
            final_protein_List.append(make_Protein)
            i=i+(3*len(make_DnaToRna))
        else:
            i=i+1
            count=count+1
        
    return final_protein_List


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    unique_Proteins=[]
    for each in proteinList1:
        #print(each)
        if each in proteinList2 and each not in unique_Proteins:
            unique_Proteins.append(each)

    return unique_Proteins


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    aminoAcids_List=[]
    for each in proteinList:
        for i in each:
            aminoAcids_List.append(i)
    #print(aminoAcids_List)

    return aminoAcids_List


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    aminoAcid_Dict={}
    #print(aaList)
    for each in aaList:
        if each not in aminoAcid_Dict:
            aminoAcid_Dict[each]= aaList.count(each)
    return aminoAcid_Dict


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    three_Element_List=[]
    #print(proteinList1)
    proteins_List1=combineProteins(proteinList1)
    aminoAcid_Dict1=aminoAcidDictionary(proteins_List1)

    proteins_List2=combineProteins(proteinList2)
    aminoAcid_Dict2=aminoAcidDictionary(proteins_List2)
    #print(cutoff)

    list1_Difference=list(set(proteins_List2) - set(proteins_List1))
    list2_Difference=list(set(proteins_List1) - set(proteins_List2))

    #print(list2_Difference)
    for each1 in list1_Difference:
        aminoAcid_Dict1[each1]=0
    for each2 in list2_Difference:
        aminoAcid_Dict2[each2]=0
        #print(aminoAcid_Dict2)
    
    #frequencies
    length1=len(proteins_List1)
    aminoAcid_Freq1={}
    for each_acid1 in aminoAcid_Dict1:
        aminoAcid_Freq1[each_acid1]=aminoAcid_Dict1[each_acid1]/length1

    length2=len(proteins_List2)
    aminoAcid_Freq2={}
    for each_acid2 in aminoAcid_Dict2:
        aminoAcid_Freq2[each_acid2]=aminoAcid_Dict2[each_acid2]/length2
    
    #result
    for key,value in aminoAcid_Dict1.items():
        innerList=[]
        if key!="Start" and key!="Stop":
            if key in aminoAcid_Dict2.keys():
                diff_Freq=abs(aminoAcid_Freq1[key]-aminoAcid_Freq2[key])
                if diff_Freq>cutoff:
                    innerList.append(key)
                    innerList.append(aminoAcid_Freq1[key])
                    innerList.append(aminoAcid_Freq2[key])
                    three_Element_List.append(innerList)

    return three_Element_List


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print(commonalities)
    #protein=""
    print("\n","These are the common proteins","\n")
    for each in commonalities:
        for i in each:
            if i!="Start" and i!="Stop":
                print(i,"\n")
   
    print("\n","These are the amino acids that occurred at the most different rates","\n")
    #print(differences)
    for each_List in differences:
        a_number1 = each_List[1]
        a_number2 = each_List[2]
        percentage1 = "{:.2%}".format(a_number1)
        percentage2 = "{:.2%}".format(a_number2)
        print(each_List[0],":",percentage1,"in Seq1,",percentage2,"in Seq2")
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    #print(proteinList1)
    amino_List1=combineProteins(proteinList1)
    list1=aminoAcidDictionary(amino_List1)
    amino_List2=combineProteins(proteinList2)
    list2=aminoAcidDictionary(amino_List2)

    #print(list1)
    label1=list(list1.keys())
    label2=list(list2.keys())

    unique_List=list(set(label2)-set(label1))
    sorted_Aminoacid_List=sorted(label1+unique_List)


    return sorted_Aminoacid_List


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    frequency_List=[]
    L=combineProteins(proteinList)
    #print(L)
    D=aminoAcidDictionary(L)
    for i in labels:
        if i in L:
            frequency=D[i]/len(L)
            frequency_List.append(frequency)
        else:
            frequency_List.append(0)

    return frequency_List


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    
    w= 0.4
    plt.bar(xLabels,freqList1,width=-w,align="edge" ,label=label1, edgecolor= edgeList)
    plt.bar(xLabels,freqList2,width=w,align="edge", label=label2, edgecolor= edgeList)
    plt.xticks(rotation="vertical")
    plt.legend()
    plt.title("comparing the two gene frequency lists")

    plt.show()
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    edge_List=[]
    diff_List=[]
    for each in biggestDiffs:
        diff_List.append(each[0])
    for each_aa in labels:
        if each_aa in diff_List:
            edge_List.append("black")
        else:
            edge_List.append("white")

    return edge_List


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # test.testReadFile()
    # test.testDnaToRna()
    # test.testMakeCodonDictionary()
    # test.testGenerateProtein()
    # test.testSynthesizeProteins()
    # test.testCommonProteins()
    # test.testCombineProteins()
    # test.testAminoAcidDictionary()
    # test.testFindAminoAcidDifferences()
    
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()


    ## Uncomment these for Week 2 ##
    
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()
    

    ## Uncomment these for Week 3 ##
    
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    

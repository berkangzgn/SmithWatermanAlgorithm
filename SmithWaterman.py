#Libraries
import numpy
import argparse
import os
import re
import sys
import unittest
import time
from pandas import *

dosya = open('cikti.txt' , 'w')

#Initialization
seqA = 'TTAGCTGATCTTAC'
seqB = 'TTAGGCTATCGA'

#Dosyadan açma işlemleri aşağıdaki gibi yapılabilir.
# seqA = open("FASTA_1.txt", "r")
# seqB = open("FASTA_2.txt", "r")

match = 3
mismatch = -1
alphabet = 'ACGT'

def main():
    #First we display the Substitution matrix, the scoring matrix
    print()
    print('Alfabe için atama matrisi {0} :'.format(alphabet))
    Substi = create_Substi(alphabet, match, mismatch)
    print()
    print(DataFrame(Substi))

    rows = len(seqA) + 1
    cols = len(seqB) + 1

    score_matrix, start_pos, max_score, path_matrix = create_score_matrix(rows, cols)

    print()
    print('Puanlama Matrisi : \n')
    print_matrix(score_matrix)
    print()
    print('En Yuksek Puan = ', max_score)
    print()
    print('start_pos=',start_pos)
    print()
    print('0a giden yol:')
    time.sleep(0.5)

    # Traceback : it also displays the path
    seqA_aligned, seqB_aligned = traceback(path_matrix, start_pos, score_matrix)
    assert len(seqA_aligned) == len(seqB_aligned), 'Hizalanmis dizeler ayni boyutta degil ! '

    # Printing the alignements and a few figures
    Graphical_display(seqA_aligned, seqB_aligned)
    print('Hizalama puani = ' ,max_score)
    dosya.write('\n\nHizalama puani = {0}'.format(max_score))
    return(0)

def traceback(path_matrix, start_pos, score_matrix):
    x,y = start_pos
    aligned_seqA = []
    aligned_seqB = []

    while path_matrix[x][y] != [0, 'NULL'] :
        d, direction = path_matrix[x][y][0], path_matrix[x][y][1]
        if direction == 'Capraz' :
            assert d ==1, 'path_matrix yanlis olusturulmus !'
            aligned_seqA.append(seqA[x - 1])
            aligned_seqB.append(seqB[y - 1])
            x -= 1
            y -= 1
            print('Capraz',score_matrix[x][y])
        elif direction == 'Yukari' :
            for c in range(d):
                aligned_seqA.append(seqA[x - 1])
                aligned_seqB.append('-')
                x -= 1
                print('Yukari',score_matrix[x][y])
        elif direction == 'Sol' :
            for c in range(d):
                aligned_seqA.append('-')
                aligned_seqB.append(seqB[y - 1])
                y -= 1
                print('Sol',score_matrix[x][y])
    print('Traceback su konumda 0 a ulasti :',(x,y))

    return ''.join(reversed(aligned_seqA)), ''.join(reversed(aligned_seqB))

def create_score_matrix(rows, cols):
    score_matrix = [[0 for col in range(cols)] for row in range(rows)]
    path_matrix = [[[0 , 'NULL'] for col in range(cols)] for row in range(rows)]

    max_score = 0
    max_pos   = None
    for i in range(1, rows):
        for j in range(1, cols):
            score, antecedent = calc_score(score_matrix, i, j)
            if score > max_score:
                max_score = score
                max_pos   = (i, j)

            score_matrix[i][j],path_matrix[i][j] = score, antecedent

    assert max_pos is not None, 'Maksimum bulunamadi :( '

    return score_matrix, max_pos, max_score, path_matrix

def calc_score(score_matrix, x, y):
    similarity = Substitution_score(Substi, x, y)

    same_row = [(score_matrix[x][y-l]-gap_penalty(l)) for l in range(1,x+1)]
    same_col = [(score_matrix[x-k][y]-gap_penalty(k)) for k in range(1,x+1)]

    up_score = max(same_col)
    left_score = max(same_row)

    diag_score = score_matrix[x - 1][y - 1] + similarity
    pos_max_up = first_pos_max(same_col)
    pos_max_left = first_pos_max(same_row)

    score =  max(0, diag_score, up_score, left_score)

    if score == diag_score :
        antecedent = [1, 'Capraz']
        return score, antecedent
    elif score == up_score :
        antecedent = [pos_max_up + 1, 'Yukari']
        return score, antecedent
    elif score == left_score :
        antecedent = [pos_max_left + 1, 'Sol']
        return score, antecedent
    else :
        return score, [0, 'NULL']

def alignment_string(aligned_seqA, aligned_seqB):

    idents, gaps, mismatches = 0, 0, 0
    alignment_string = []
    for base1, base2 in zip(aligned_seqA, aligned_seqB): #loop that runs over both sequences simultaneously
        if base1 == base2:
            alignment_string.append('|')
            idents += 1
        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1
        else:
            alignment_string.append(':')
            mismatches += 1

    return ''.join(alignment_string), idents, gaps, mismatches

def create_Substi(alphabet, match, mismatch):
    global Substi
    Substi = [['NULL' for col in range(len(alphabet))] for row in range(len(alphabet))]

    for i in range(len(alphabet)):
        for j in range(len(alphabet)):
            if alphabet[i] == alphabet[j] :
                Substi[i][j] = match
            elif alphabet[i] != alphabet[j]:
                Substi[i][j] = mismatch

    return Substi

def Substitution_score(Substi, x, y):
    a_i = alphabet.index(seqA[x-1])
    b_j = alphabet.index(seqB[y-1])

    return Substi[a_i][b_j]

def Graphical_display(seqA_aligned, seqB_aligned):
    alignment_str, idents, gaps, mismatches = alignment_string(seqA_aligned, seqB_aligned)
    alength = len(seqA_aligned)
    time.sleep(0.2)
    print()
    print('Eslesme orani = {0}/{1} ({2:.1%})'.format(idents,
          alength, idents / alength))
    dosya.write("Eslesme orani = {0}/{1} ({2:.1%})".format(idents,alength, idents / alength))
    print('Eslesmeme orani = {0}/{1} ({2:.1%})'.format(mismatches, alength, mismatches / alength))
    dosya.write('\nEslesmeme orani = {0}/{1} ({2:.1%})'.format(mismatches, alength, mismatches / alength))
    print('Bosluk orani = {0}/{1} ({2:.1%})'.format(gaps, alength, gaps / alength))
    dosya.write('\nBosluk orani = {0}/{1} ({2:.1%})'.format(gaps, alength, gaps / alength))
    print()

    for i in range(0, alength, 60):
        seqA_slice = seqA_aligned[i:i+60]
        print('Gen_1 uzunluk:  {0:<4}   '.format( i + len(seqA_slice) , seqA_slice ))
        dosya.write('\n\nGen_1 uzunluk:  {0:<4}   '.format( i + len(seqA_slice) , seqA_slice ))

        seqB_slice = seqB_aligned[i:i+60]
        print('Gen_2 uzunluk:  {0:<4}  '.format(i + len(seqB_slice) , seqB_slice ))
        dosya.write('\nGen_2 uzunluk:  {0:<4}   '.format(i + len(seqB_slice) , seqB_slice ))

        dosya.write('\n\nHizalama Sonucu:')
        dosya.write('\n{0}'.format(seqA_slice))
        dosya.write('\n{0}'.format(alignment_str[i:i+60]))
        dosya.write('\n{0}'.format(seqA_slice))

        print('')
        print('Hizalama Sonucu:')
        print('{0}'.format(seqA_slice))
        print('{0}'.format(alignment_str[i:i+60]))
        print('{0}'.format(seqB_slice))
        print()

def print_matrix(matrix): # graphical display of matrix
    print('\n'.join([''.join(['     {:4}'.format(item) for item in row])for row in matrix]))

def first_pos_max(list):
    maxi = max(list)
    return [i for i, j in enumerate(list) if j == maxi][0]

def gap_penalty(k) :
    return 3*k

if __name__ == '__main__':
    sys.exit(main())

dosya.close()

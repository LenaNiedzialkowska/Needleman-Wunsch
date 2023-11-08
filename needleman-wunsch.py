#!/usr/bin/env python3
import numpy as np
from Bio import SeqIO


# Pobranie z pliku .fasta sekwencji
seq_1 = True
input_file = input("Podaj plik fasta z sekwencjami: ")
fasta_sequences = SeqIO.parse(open(input_file),'fasta')
for fasta in fasta_sequences:
    if seq_1 == True:
        name_1, sequence_1 = fasta.id, str(fasta.seq)
        seq_1 = False
    else:
        name_2, sequence_2 = fasta.id, str(fasta.seq)

# print(sequence_1, sequence_2)

# Stworzenie macierzy o wymiarach (n+1)x(n+1)
matrix = np.zeros((len(sequence_1) + 1, len(sequence_2) + 1))
check_matrix_score = np.zeros((len(sequence_1), len(sequence_2)))

match_score = float(input("Podaj match score: "))
mismatch_score = float(input("Podaj mismatch score: "))
gap_score = float(input("Podaj gap score: "))

# Wypełnianie macierzy sprawdzającej wartościami odpowiadającymi match_score, mismatch_score i gap_score
for i in range(len(sequence_1)):
    for j in range(len(sequence_2)):
        if sequence_1[i] == sequence_2[j]:
            check_matrix_score[i][j] = match_score
        else:
            check_matrix_score[i][j] = mismatch_score


# Algorytm Needlemana-Wunscha

# Inicjalizacja macierzy
for i in range(len(sequence_1) + 1):
    matrix[i][0] = i * gap_score
for j in range(len(sequence_2) + 1):
    matrix[0][j] = j * gap_score

# Uzupełnianie macierzy zgodnie z algorytmem
for i in range(1, len(sequence_1) + 1):
    for j in range(1, len(sequence_2) + 1):
        matrix[i][j] = max(matrix[i - 1][j - 1] + check_matrix_score[i - 1][j - 1],
                           matrix[i - 1][j] + gap_score,
                           matrix[i][j - 1] + gap_score)


# Dopasowanie sekwencji
aligned_1 = ""
aligned_2 = ""

len_1 = len(sequence_1)
len_2 = len(sequence_2)

while len_1 > 0 and len_2 > 0:
    if len_1 > 0 and len_2 > 0 and matrix[len_1][len_2] == matrix[len_1 - 1][len_2 - 1] + check_matrix_score[len_1 - 1][len_2 - 1]:
        aligned_1 = sequence_1[len_1 - 1] + aligned_1
        aligned_2 = sequence_2[len_2 - 1] + aligned_2

        len_1 -= 1
        len_2 -= 1

    elif len_1 > 0 and matrix[len_1][len_2] == matrix[len_1 - 1][len_2] + gap_score:
        aligned_1 = sequence_1[len_1 - 1] + aligned_1
        aligned_2 = "-" + aligned_2

        len_1 -= 1

    else:
        aligned_1 = "-" + aligned_1
        aligned_2 = sequence_2[len_2 - 1] + aligned_2

        len_2 -= 1


output_file = open("output.txt", "w")
output_file.write(aligned_1)
output_file.write('\n')
output_file.write(aligned_2)
output_file.close()

# print(aligned_1)
# print(aligned_2)
# print('Score: ', matrix[len(sequence_1)][len(sequence_2)])

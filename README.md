# Implementacja algorytmu Needlemana-Wunscha w pythonie
Program przyjmuje na wejściu plik w formacie fasta (np. sample.fasta) oraz 3 argumenty: match_score, mismatch_score i gap_score.

W celu uruchomienia programu należy podać komendę:
python needleman-wunsch.py 

<ścieżka do pliku>
<match_score>
<mismatch_score>
<gap_score>

Dla argumentów match_score = 1, mismatch_score = -1 i gap_score = -2 oraz sekwencji: 

>DOJHLOP01DWJWJ length=100 xy=1483_3377 region=1 run=R_2005_09_08_15_35_38_
TCGATTCTATGGAGGGATGCTGGCAAGGCTCCGGAAGCAGCATCAGCAATTAAAAAATTA
CTGGACCTGATCTTATGAAGTTAGGATTGTTGACGAGGTA
>DOJHLOP01CAI84 length=98 xy=0822_3942 region=1 run=R_2005_09_08_15_35_38_
AGGCGTCGCAGACAGGTTACTTATGTTTGAACATAGTGTTTACACAGTTGCAAGCCCTGA
AGCTTGTGCTTCGATTCTATGGAGGGATGCTGGCAAGG

Wynikiem będzie:

TCGATTCTATGGAGGGATGCTGGCAAGGCTCCGGA-AGCAGCAT-CAGCAATTAAAAAATTACTGGA-CCTGATCTT--A-TGAAGTTAGGATTGTTGACGAGGTA

AGGCGTCGCAGACAGGTTACT--TATGTTTGAACATAG-TGTTTACA-CAGTT--GCAAGCCCTGAAGCTTGTGCTTCGATTCTATGGAGGGATGCTGGCAAGG--

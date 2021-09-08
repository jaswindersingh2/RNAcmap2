#Name: PKB246
#Length:  34 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GACCGTGCACCATGAGCACGAATCCTAAACCTCA
.[.............(...]....)...[(..])
EEEEEEEEEEEEEEESHHHHHHHHSXXXXSHHHS
NKNNNNNNNNNNNNNNNNNKNNNNNNNNKNNNKN
S1 16..16 "G" 25..25 "C"
S2 30..30 "C" 34..34 "A"
H1 17..24 "CACGAATC" (16,25) G:C PK{1}
H2 31..33 "CTC" (30,34) C:A PK{2}
X1 26..29 "TAAA" (25,16) C:G (30,34) C:A PK{2}
E1 1..15 "GACCGTGCACCATGA" PK{1}
PK1 1bp 2..2 20..20 E1 1..15 H1 17..24
PK2 1bp 29..29 33..33 X1 26..29 H2 31..33
PK1.1 2 A 20 G
PK2.1 29 A 33 C
NCBP1 30 C 34 A S2
NCBP2 2 A 20 G PK1.1
NCBP3 29 A 33 C PK2.1
segment1 1bp 16..16 G 25..25 C
segment2 1bp 30..30 C 34..34 A

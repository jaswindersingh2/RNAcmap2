#Name: PKB246
#Length:  34 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GACCGTGCACCATGAGCACGAATCCTAAACCTCA
[......(.{].......(..}).)...[(..])
EEEEEEESIIIIIIIIIISHHHSISXXXXSHHHS
KNNNNNNNNKKNNNNNNNNNNKNNNNNNKNNNKN
S1 8..8 "C" 25..25 "C"
S2 19..19 "C" 23..23 "T"
S3 30..30 "C" 34..34 "A"
H1 20..22 "GAA" (19,23) C:T PK{2}
H2 31..33 "CTC" (30,34) C:A PK{3}
I1.1 9..18 "ACCATGAGCA" (8,25) C:C PK{1,2}
I1.2 24..24 "C" (23,19) T:C 
X1 26..29 "TAAA" (25,8) C:C (30,34) C:A PK{3}
E1 1..7 "GACCGTG" PK{1}
PK1 1bp 1..1 11..11 E1 1..7 I1.1 9..18
PK2 1bp 10..10 22..22 I1.1 9..18 H1 20..22
PK3 1bp 29..29 33..33 X1 26..29 H2 31..33
PK1.1 1 G 11 C
PK2.1 10 C 22 A
PK3.1 29 A 33 C
NCBP1 8 C 25 C S1
NCBP2 19 C 23 T S2
NCBP3 30 C 34 A S3
NCBP4 10 C 22 A PK2.1
NCBP5 29 A 33 C PK3.1
segment1 2bp 8..19 CACCATGAGCAC 23..25 TCC
segment2 1bp 30..30 C 34..34 A

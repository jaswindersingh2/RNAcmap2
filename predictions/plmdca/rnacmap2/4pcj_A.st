#Name: 4pcj_A
#Length:  35 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
CUGCUGGCUAAGGCAUGAAAGUGCUAUGCCUGCUG
[(.....)......{.(...]}..(....))....
ESHHHHHSXXXXXXXXSBBBBBBBSHHHHSSEEEE
KNNNNNNNNNNNNNKNNNNNKKNNNNNNNNNNNNN
S1 2..2 "U" 8..8 "C"
S2 17..17 "G" 31..31 "U"
S3 25..25 "U" 30..30 "C"
H1 3..7 "GCUGG" (2,8) U:C 
H2 26..29 "AUGC" (25,30) U:C 
B1 18..24 "AAAGUGC" (17,31) G:U (25,30) U:C PK{1,2}
X1 9..16 "UAAGGCAU" (8,2) C:U (17,31) G:U PK{2}
E1 1..1 "C" PK{1}
E2 32..35 "GCUG" 
PK1 1bp 1..1 21..21 E1 1..1 B1 18..24
PK2 1bp 15..15 22..22 X1 9..16 B1 18..24
PK1.1 1 C 21 G
PK2.1 15 A 22 U
NCBP1 2 U 8 C S1
NCBP2 25 U 30 C S3
segment1 1bp 2..2 U 8..8 C
segment2 2bp 17..25 GAAAGUGCU 30..31 CU

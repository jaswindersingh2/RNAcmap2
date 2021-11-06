#Name: 1s03_A
#Length: ,47,
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGACGAUGGCGAAACUGCAUGAGGCAAUUCAUGCAAGUCCCUCGUCC
....[..(..........[..(]........)..........]...)
EEEEEEESIIIIIIIIIIIIISHHHHHHHHHSIIIIIIIIIIIIIIS
NNNNKNNNNNNNNNNNNNKNNNKNNNNNNNNNNNNNNNNNNNKNNNN
S1 8..8 "G" 47..47 "C"
S2 22..22 "A" 32..32 "U"
H1 23..31 "GGCAAUUCA" (22,32) A:U PK{2}
I1.1 9..21 "GCGAAACUGCAUG" (8,47) G:C PK{2}
I1.2 33..46 "GCAAGUCCCUCGUC" (32,22) U:A PK{1}
E1 1..7 "GGACGAU" PK{1}
PK1 1bp 5..5 43..43 E1 1..7 I1.2 33..46
PK2 1bp 19..19 23..23 I1.1 9..21 H1 23..31
PK1.1 5 G 43 C
PK2.1 19 A 23 G
NCBP1 19 A 23 G PK2.1
segment1 2bp 8..22 GGCGAAACUGCAUGA 32..47 UGCAAGUCCCUCGUCC

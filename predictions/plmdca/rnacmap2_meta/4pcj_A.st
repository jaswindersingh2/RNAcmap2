#Name: 4pcj_A
#Length: ,35,
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
CUGCUGGCUAAGGCAUGAAAGUGCUAUGCCUGCUG
[..(...{..(..((.....])).}).....)...
EEESIIIIIISIISSHHHHHHSSIISIIIIISEEE
KNNNNNNKNNNNNNNNNNNNKNNNKNNNNNNNNNN
S1 4..4 "C" 32..32 "G"
S2 11..11 "A" 26..26 "A"
S3 14..15 "CA" 22..23 "UG"
H1 16..21 "UGAAAG" (15,22) A:U PK{1}
I1.1 5..10 "UGGCUA" (4,32) C:G PK{2}
I1.2 27..31 "UGCCU" (26,11) A:A 
I2.1 12..13 "GG" (11,26) A:A 
I2.2 24..25 "CU" (23,14) G:C PK{2}
E1 1..3 "CUG" PK{1}
E2 33..35 "CUG" 
PK1 1bp 1..1 21..21 E1 1..3 H1 16..21
PK2 1bp 8..8 25..25 I1.1 5..10 I2.2 24..25
PK1.1 1 C 21 G
PK2.1 8 C 25 U
NCBP1 11 A 26 A S2
NCBP2 8 C 25 U PK2.1
segment1 4bp 4..15 CUGGCUAAGGCA 22..32 UGCUAUGCCUG

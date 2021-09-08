#Name: PKB246
#Length:  34 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GACCGTGCACCATGAGCACGAATCCTAAACCTCA
([..).....(...)(...]....)...[(..])
SHHHSMMMMMSHHHSSHHHHHHHHSMMMMSHHHS
NKNNNNNNNNNNNNNNNNNKNNNNNNNNKNNNKN
S1 1..1 "G" 5..5 "G"
S2 11..11 "C" 15..15 "A"
S3 16..16 "G" 25..25 "C"
S4 30..30 "C" 34..34 "A"
H1 2..4 "ACC" (1,5) G:G PK{1}
H2 12..14 "ATG" (11,15) C:A 
H3 17..24 "CACGAATC" (16,25) G:C PK{1}
H4 31..33 "CTC" (30,34) C:A PK{2}
M1.1 6..10 "TGCAC" (5,1) G:G (11,15) C:A 
M1.2 16..15 "" (15,11) A:C (16,25) G:C 
M1.3 26..29 "TAAA" (25,16) C:G (30,34) C:A PK{2}
PK1 1bp 2..2 20..20 H1 2..4 H3 17..24
PK2 1bp 29..29 33..33 M1.3 26..29 H4 31..33
PK1.1 2 A 20 G
PK2.1 29 A 33 C
NCBP1 1 G 5 G S1
NCBP2 30 C 34 A S4
NCBP3 11 C 15 A S2
NCBP4 2 A 20 G PK1.1
NCBP5 29 A 33 C PK2.1
segment1 1bp 1..1 G 5..5 G
segment2 1bp 11..11 C 15..15 A
segment3 1bp 16..16 G 25..25 C
segment4 1bp 30..30 C 34..34 A

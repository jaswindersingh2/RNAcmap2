#Name: 4x4v_B
#Length:  37 
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGCGCUGCGGGGUUCGAGUCCCCGCAGUGUUGCCACC
.....[.........(.]..)....(....)(....)
EEEEEEEEEEEEEEESHHHHSXXXXSHHHHSSHHHHS
NNNNNKNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNN
S1 16..16 "G" 21..21 "C"
S2 26..26 "A" 31..31 "U"
S3 32..32 "G" 37..37 "C"
H1 17..20 "AGUC" (16,21) G:C PK{1}
H2 27..30 "GUGU" (26,31) A:U 
H3 33..36 "CCAC" (32,37) G:C 
X1 22..25 "CCGC" (21,16) C:G (26,31) A:U 
E1 1..15 "GGCGCUGCGGGGUUC" PK{1}
PK1 1bp 6..6 18..18 E1 1..15 H1 17..20
PK1.1 6 U 18 G
segment1 1bp 16..16 G 21..21 C
segment2 1bp 26..26 A 31..31 U
segment3 1bp 32..32 G 37..37 C

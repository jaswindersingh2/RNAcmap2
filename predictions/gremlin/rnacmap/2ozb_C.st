#Name: 2ozb_C
#Length:  33 
#PageNumber: 3
AUCGUAGCCAAUGAGGUUUAUCCGAGGCGCGAU
([{{)..(.......(.....)]...)(.}}.)
SHHHSMMSIIIIIIISHHHHHSIIIISSHHHHS
NKKKNNNNNNNNNNNNNNNNNNKNNNNNNKKNN
S1 1..1 "A" 5..5 "U"
S2 8..8 "C" 27..27 "G"
S3 16..16 "G" 22..22 "C"
S4 28..28 "C" 33..33 "U"
H1 2..4 "UCG" (1,5) A:U PK{1,2}
H2 17..21 "UUUAU" (16,22) G:C 
H3 29..32 "GCGA" (28,33) C:U PK{2}
I1.1 9..15 "CAAUGAG" (8,27) C:G 
I1.2 23..26 "CGAG" (22,16) C:G PK{1}
M1.1 6..7 "AG" (5,1) U:A (8,27) C:G 
M1.2 28..27 "" (27,8) G:C (28,33) C:U 
PK1 1bp 2..2 23..23 H1 2..4 I1.2 23..26
PK2 2bp 3..4 30..31 H1 2..4 H3 29..32
PK1.1 2 U 23 C
PK2.1 3 C 31 G
PK2.2 4 G 30 C
NCBP1 28 C 33 U S4
NCBP2 2 U 23 C PK1.1
segment1 1bp 1..1 A 5..5 U
segment2 2bp 8..16 CCAAUGAGG 22..27 CCGAGG
segment3 1bp 28..28 C 33..33 U

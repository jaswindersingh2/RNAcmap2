#Name: 5x2h_B
#Length:  93 
#PageNumber: 2
GGAAAUUAGGUGCGCUUGGCGUUUUAGUCCCUGAAAAGGGACUAAAAUAAAGAGUUUGCGGGACUCUGCGGGGUUACAAUCCCCUAAAACCGC
[(.....)............(((((((((((.]...).))))))))))..((..(((...)).)).)..(..(......[)..)]...(...)
ESHHHHHSXXXXXXXXXXXXSSSSSSSSSSSHHHHHSBSSSSSSSSSSXXSSBBSSSHHHSSBSSBSXXSIISHHHHHHHSIISXXXXSHHHS
KNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNKNNNNNNNN
S1 2..2 "G" 8..8 "A"
S2 21..30 "GUUUUAGUCC" 39..48 "GGACUAAAAU"
S3 31..31 "C" 37..37 "A"
S4 51..51 "A" 67..67 "U"
S5 52..52 "G" 65..65 "U"
S6 55..55 "U" 64..64 "C"
S7 56..57 "UU" 61..62 "GG"
S8 70..70 "G" 84..84 "C"
S9 73..73 "G" 81..81 "C"
S10 89..89 "A" 93..93 "C"
H1 3..7 "AAAUU" (2,8) G:A 
H2 32..36 "UGAAA" (31,37) C:A PK{1}
H3 58..60 "GCG" (57,61) U:G 
H4 74..80 "UUACAAU" (73,81) G:C PK{2}
H5 90..92 "CCG" (89,93) A:C 
B1 38..38 "G" (37,31) A:C (39,30) G:C 
B2 53..54 "AG" (52,65) G:U (55,64) U:C 
B3 63..63 "A" (62,56) G:U (64,55) C:U 
B4 66..66 "C" (65,52) U:G (67,51) U:A 
I1.1 71..72 "GG" (70,84) G:C 
I1.2 82..83 "CC" (81,73) C:G 
X1 9..20 "GGUGCGCUUGGC" (8,2) A:G (21,48) G:U 
X2 49..50 "AA" (48,21) U:G (51,67) A:U 
X3 68..69 "GC" (67,51) U:A (70,84) G:C 
X4 85..88 "UAAA" (84,70) C:G (89,93) A:C PK{2}
E1 1..1 "G" PK{1}
PK1 1bp 1..1 33..33 E1 1..1 H2 32..36
PK2 1bp 80..80 85..85 H4 74..80 X4 85..88
PK1.1 1 G 33 G
PK2.1 80 U 85 U
NCBP1 55 U 64 C S6
NCBP2 89 A 93 C S10
NCBP3 31 C 37 A S3
NCBP4 2 G 8 A S1
NCBP5 1 G 33 G PK1.1
NCBP6 80 U 85 U PK2.1
segment1 1bp 2..2 G 8..8 A
segment2 11bp 21..31 GUUUUAGUCCC 37..48 AGGGACUAAAAU
segment3 5bp 51..57 AGAGUUU 61..67 GGACUCU
segment4 2bp 70..73 GGGG 81..84 CCCC
segment5 1bp 89..89 A 93..93 C

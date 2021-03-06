#Name: 5b2o_B
#Length:  94 
#PageNumber: 3
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGGAAAUUAGGUGCGCUGGGGGUUUCAGUUGCGCCGAAAGGCGCUCUGUAAUCAUUUAAAAGUAUUUUGAACGGACCUCUGUUUGACACGUCUG
.....(..................[........................(.......{..].)....(......[..).....).......}.]
EEEEESMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMSHHHHHHHHHHHHSMMMMSHHHHHHHHHSMMMMMSEEEEEEEEEE
NNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNKNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNKNK
S1 6..6 "A" 84..84 "U"
S2 50..50 "A" 63..63 "U"
S3 68..68 "U" 78..78 "U"
H1 51..62 "AUCAUUUAAAAG" (50,63) A:U PK{1,2}
H2 69..77 "GAACGGACC" (68,78) U:U PK{3}
M1.1 7..49 "UUAGGUGCGCUGGGGGUUUCAGUUGCGCCGAAAGGCGCUCUGU" (6,84) A:U (50,63) A:U PK{1}
M1.2 64..67 "AUUU" (63,50) U:A (68,78) U:U 
M1.3 79..83 "CUGUU" (78,68) U:U (84,6) U:A 
E1 1..5 "GGGAA" 
E2 85..94 "GACACGUCUG" PK{2,3}
PK1 1bp 25..25 61..61 M1.1 7..49 H1 51..62
PK2 1bp 58..58 92..92 H1 51..62 E2 85..94
PK3 1bp 75..75 94..94 H2 69..77 E2 85..94
PK1.1 25 U 61 A
PK2.1 58 A 92 C
PK3.1 75 A 94 G
NCBP1 68 U 78 U S3
NCBP2 58 A 92 C PK2.1
NCBP3 75 A 94 G PK3.1
segment1 1bp 6..6 A 84..84 U
segment2 1bp 50..50 A 63..63 U
segment3 1bp 68..68 U 78..78 U

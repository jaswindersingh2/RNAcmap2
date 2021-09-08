#Name: 6tfe_A
#Length: ,52,
#PageNumber: 2
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
#Warning: Structure contains linked PK-segments of same sizes 1 and 1. Using PK brackets for the more 5' segment
GGCUUCAACAACCCCGUAGGUUGGGCCGAAAGGCAGCGAAUCUACUGGAGCC
...[......(.....[(.......((....))...].....)]...)....
EEEEEEEEEESIIIIIISIIIIIIISSHHHHSSIIIIIIIIISIIIISEEEE
NNNKNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNKNNNNNNKNNNNNNNN
S1 11..11 "A" 48..48 "G"
S2 18..18 "A" 43..43 "U"
S3 26..27 "CC" 32..33 "GG"
H1 28..31 "GAAA" (27,32) C:G 
I1.1 12..17 "CCCCGU" (11,48) A:G PK{2}
I1.2 44..47 "ACUG" (43,18) U:A PK{1}
I2.1 19..25 "GGUUGGG" (18,43) A:U 
I2.2 34..42 "CAGCGAAUC" (33,26) G:C PK{2}
E1 1..10 "GGCUUCAACA" PK{1}
E2 49..52 "AGCC" 
PK1 1bp 4..4 44..44 E1 1..10 I1.2 44..47
PK2 1bp 17..17 37..37 I1.1 12..17 I2.2 34..42
PK1.1 4 U 44 A
PK2.1 17 U 37 C
NCBP1 11 A 48 G S1
NCBP2 17 U 37 C PK2.1
segment1 4bp 11..27 ACCCCGUAGGUUGGGCC 32..48 GGCAGCGAAUCUACUGG

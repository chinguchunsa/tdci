V33 :0x4 get_ham0
15 mod_getham0.f90 S624 0
11/27/2019  13:06:56
use units public 0 indirect
use setup_variables public 0 indirect
use route_control public 0 indirect
use global_variables public 0 direct
enduse
S 624 24 0 0 0 8 1 0 5015 10005 0 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 get_ham0
S 1415 6 4 0 0 7 1416 624 11763 4 0 A 0 0 0 0 B 0 7 0 0 0 0 0 0 0 0 0 0 1419 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 double_startaa
S 1416 6 4 0 0 7 1417 624 11778 4 0 A 0 0 0 0 B 0 7 0 0 0 8 0 0 0 0 0 0 1419 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 double_startbb
S 1417 6 4 0 0 7 1418 624 11793 4 0 A 0 0 0 0 B 0 7 0 0 0 16 0 0 0 0 0 0 1419 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 double_startab
S 1418 6 4 0 0 7 1 624 11808 4 0 A 0 0 0 0 B 0 7 0 0 0 24 0 0 0 0 0 0 1419 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 double_startba
S 1419 11 0 0 0 8 1401 624 11823 40800000 805000 A 0 0 0 0 B 0 19 0 0 0 32 0 0 1415 1418 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _get_ham0$2
S 1420 23 5 0 0 0 1421 624 11835 0 0 A 0 0 0 0 B 0 133 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_cis
S 1421 14 5 0 0 0 1 1420 11835 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 23 0 0 0 0 0 0 0 0 0 0 0 0 0 23 0 624 0 0 0 0 get_cis
F 1421 0
S 1422 23 5 0 0 0 1423 624 11843 0 0 A 0 0 0 0 B 0 381 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_soc_cis
S 1423 14 5 0 0 0 1 1422 11843 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 24 0 0 0 0 0 0 0 0 0 0 0 0 0 137 0 624 0 0 0 0 get_soc_cis
F 1423 0
S 1424 23 5 0 0 0 1425 624 11855 0 0 A 0 0 0 0 B 0 478 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_soc_ao2mo
S 1425 14 5 0 0 0 1 1424 11855 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 25 0 0 0 0 0 0 0 0 0 0 0 0 0 385 0 624 0 0 0 0 get_soc_ao2mo
F 1425 0
S 1426 23 5 0 0 0 1427 624 11869 0 0 A 0 0 0 0 B 0 541 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_cis_index
S 1427 14 5 0 0 0 1 1426 11869 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 482 0 624 0 0 0 0 get_cis_index
F 1427 0
S 1428 23 5 0 0 0 1429 624 11883 0 0 A 0 0 0 0 B 0 720 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_ip_cisd
S 1429 14 5 0 0 0 1 1428 11883 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 27 0 0 0 0 0 0 0 0 0 0 0 0 0 545 0 624 0 0 0 0 get_ip_cisd
F 1429 0
S 1430 23 5 0 0 0 1431 624 11895 0 0 A 0 0 0 0 B 0 806 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 get_ip_index
S 1431 14 5 0 0 0 1 1430 11895 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 28 0 0 0 0 0 0 0 0 0 0 0 0 0 724 0 624 0 0 0 0 get_ip_index
F 1431 0
S 1432 23 5 0 0 0 1433 624 11908 0 0 A 0 0 0 0 B 0 857 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 diagonalize
S 1433 14 5 0 0 0 1 1432 11908 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 810 0 624 0 0 0 0 diagonalize
F 1433 0
S 1434 23 5 0 0 0 1435 624 11920 0 0 A 0 0 0 0 B 0 910 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 zdiagonalize
S 1435 14 5 0 0 0 1 1434 11920 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 30 0 0 0 0 0 0 0 0 0 0 0 0 0 861 0 624 0 0 0 0 zdiagonalize
F 1435 0
S 1436 23 5 0 0 0 1439 624 11933 0 0 A 0 0 0 0 B 0 955 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 writeme_ham0
S 1437 1 3 1 0 28 1 1436 11946 4 43000 A 0 0 0 0 B 0 955 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 myroutine
S 1438 1 3 1 0 28 1 1436 11956 4 43000 A 0 0 0 0 B 0 955 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 option
S 1439 14 5 0 0 0 1 1436 11933 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 31 2 0 0 0 0 0 0 0 0 0 0 0 0 914 0 624 0 0 0 0 writeme_ham0
F 1439 2 1437 1438
S 1440 23 5 0 0 0 1444 624 11963 0 0 A 0 0 0 0 B 0 992 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 errors_getham0
S 1441 1 3 1 0 28 1 1440 11946 4 43000 A 0 0 0 0 B 0 992 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 myroutine
S 1442 1 3 1 0 28 1 1440 11956 4 43000 A 0 0 0 0 B 0 992 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 option
S 1443 1 3 1 0 28 1 1440 11978 80000004 43000 A 0 0 0 0 B 0 992 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mystring
S 1444 14 5 0 0 0 1 1440 11963 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 34 3 0 0 0 0 0 0 0 0 0 0 0 0 959 0 624 0 0 0 0 errors_getham0
F 1444 3 1441 1442 1443
Z
Z